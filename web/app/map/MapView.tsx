"use client";
import React, { useEffect, useMemo, useRef, useState } from "react";
import * as d3 from "d3";
import ParentView from "@/components/ParentView";
import { Box, Button, Card, Checkbox, CircularProgress, Dropdown, LinearProgress, Option, Select, Typography } from "@mui/joy";
import { useAppDispatch, useAppSelector } from "@/store/hooks";
import { View } from "@/store/features/views/types";
import TrackSelector from "./components/TrackSelector";
import CanvasHeatmap from "./components/CanvasHeatmap";
import { getViewConfig, updateViewConfig } from "@/store/features/space/spaceSlice";
import { buildTileCacheKey, TILE_SIZE_BINS, tileCache } from "../map/TileCache";
import { buildApiUrl } from "@/app/config/env";

const RESOLUTION_PRESETS = [1, 5, 10, 25, 50, 100, 200, 500, 1000];
const AUTO_TILE_CEILING = 10;
const DEFAULT_MANUAL_RESOLUTION = 100;

interface MapViewProps {
  trackFiles: any[];
  viewConfig: View;
  handleViewUpdate: (index, viewState: View) => void;
  index: number;
  dependencies: any;
}

export interface Segment {
  id: string;
  length: number;
}

const MapView = (props: MapViewProps) => {
  const [reference, setReference] = useState(props.viewConfig.config.reference);
  const [track, setTrack] = useState(props.viewConfig.config.track);
  const resolutionConfigValue = props.viewConfig.config.resolution;
  const storedManualResolution = props.viewConfig.config.manualResolution as number | undefined;
  const storedAutoResolution = props.viewConfig.config.autoResolution as number | undefined;
  const initialManualResolution =
    typeof resolutionConfigValue === "number"
      ? resolutionConfigValue
      : storedManualResolution ?? DEFAULT_MANUAL_RESOLUTION;
  const initialResolutionSetting =
    resolutionConfigValue === "auto" ? "auto" : initialManualResolution;
  const initialAutoResolution =
    typeof resolutionConfigValue === "number"
      ? resolutionConfigValue
      : storedAutoResolution ?? RESOLUTION_PRESETS[RESOLUTION_PRESETS.length - 1];
  const [resolutionSetting, setResolutionSetting] = useState<number | "auto">(initialResolutionSetting);
  const [manualResolution, setManualResolution] = useState<number>(initialManualResolution);
  const [autoResolution, setAutoResolution] = useState<number>(initialAutoResolution);
  const [availableSegments, setAvailableSegments] = useState<Segment[]>([]);
  const [toggleColourScheme, setToggleColourScheme] = useState(props.viewConfig.config.toggleColourScheme || false);
  const [normalise, setNormalise] = useState(false);
  const [negativeStrand, setNegativeStrand] = useState(props.viewConfig.config.negativeStrand || false);
  const [loading, setLoading] = useState(false);
  const [matrix, setMatrix] = useState<number[][]>([]);
  const [tileCacheVersion, setTileCacheVersion] = useState(0);
  const canvasRef = useRef<HTMLCanvasElement | null>(null);

  const heatmapRef = useRef(null);
  const space = useAppSelector((state) => state.space);
  const dispatch = useAppDispatch();

  const extractSegmentId = (value: Segment | string | null | undefined) => {
    if (!value) return null;
    if (typeof value === "string") return value;
    return value.id;
  };

  const initialSegmentAId = extractSegmentId(
    space.views?.[props.index]?.config.segmentA ?? props.viewConfig.config.segmentA,
  );
  const initialSegmentBId = extractSegmentId(
    space.views?.[props.index]?.config.segmentB ?? props.viewConfig.config.segmentB,
  );
  const [segmentAId, setSegmentAId] = useState<string | null>(initialSegmentAId);
  const [segmentBId, setSegmentBId] = useState<string | null>(initialSegmentBId);

  const segmentA = useMemo(
    () => availableSegments.find((segment) => segment.id === segmentAId) ?? null,
    [availableSegments, segmentAId],
  );
  const segmentB = useMemo(
    () => availableSegments.find((segment) => segment.id === segmentBId) ?? null,
    [availableSegments, segmentBId],
  );

  const [renderedSegmentA, setRenderedSegmentA] = useState<Segment | null>(segmentA);
  const [renderedSegmentB, setRenderedSegmentB] = useState<Segment | null>(segmentB);
  const [renderedResolution, setRenderedResolution] = useState(
    resolutionSetting === "auto" ? initialAutoResolution : initialManualResolution,
  );

  const [showTrackPicker, setShowTrackPicker] = useState(false);

  const [xLocus, setXLocus] = useState<[number, number]>([0, segmentA?.length ?? 0]);
  const [yLocus, setYLocus] = useState<[number, number]>([0, segmentB?.length ?? 0]);
  const zoomToLocusRef = useRef<((x0: number, x1: number, y0: number, y1: number) => void) | null>(null);

  useEffect(() => {
    setXLocus([0, segmentA?.length ?? 0]);
  }, [segmentA]);

  useEffect(() => {
    setYLocus([0, segmentB?.length ?? 0]);
  }, [segmentB]);

  useEffect(() => {
    const nextResolutionSetting = props.viewConfig.config.resolution;
    const nextManualResolution = props.viewConfig.config.manualResolution as number | undefined;
    const nextAutoResolution = props.viewConfig.config.autoResolution as number | undefined;
    if (nextResolutionSetting === "auto") {
      if (resolutionSetting !== "auto") {
        setResolutionSetting("auto");
      }
      if (typeof nextManualResolution === "number" && nextManualResolution !== manualResolution) {
        setManualResolution(nextManualResolution);
      }
      if (typeof nextAutoResolution === "number" && nextAutoResolution !== autoResolution) {
        setAutoResolution(nextAutoResolution);
      }
    } else if (typeof nextResolutionSetting === "number") {
      if (manualResolution !== nextResolutionSetting) {
        setManualResolution(nextResolutionSetting);
      }
      if (resolutionSetting !== nextResolutionSetting) {
        setResolutionSetting(nextResolutionSetting);
      }
      if (resolutionSetting !== "auto" && autoResolution !== nextResolutionSetting) {
        setAutoResolution(nextResolutionSetting);
      }
    }
  }, [
    props.viewConfig.config.resolution,
    props.viewConfig.config.manualResolution,
    props.viewConfig.config.autoResolution,
  ]);


  const handleTrackConfirm = (ref: string, trk: string) => {
    setReference(ref);
    setTrack(trk);
    setShowTrackPicker(false);
    props.handleViewUpdate(props.index, {
      ...props.viewConfig,
      config: {
        ...props.viewConfig.config,
        reference: ref,
        track: trk,
      },
    });
  };

  useEffect(() => {
    const fetchData = async () => {
      const queryParams = new URLSearchParams({
        uuid: space.uuid,
        file_name: track,
        genome_path: reference,
      }).toString();
        const url = new URL(buildApiUrl("/api/timge/get_segments/"));
        url.search = queryParams;
        const response = await fetch(url.toString(), {
            method: "GET",
            headers: {
            "Content-Type": "application/json",
            },
        });
      const data = await response.json();
      if (data.status === "success" && data.segments) {
        const segments: Segment[] = Object.entries(data.segments).map(([id, length]) => ({
          id,
          length: length as number,
        }));
        setAvailableSegments(segments);
      } else {
        console.error("Failed to fetch segments", data.message);
      }
    };

    if (reference && track) {
      fetchData();
    }
  }, [reference, track]);

  useEffect(() => {
    const currentConfig = props.viewConfig.config || {};
    const currentSegmentAId = extractSegmentId(currentConfig.segmentA);
    const currentSegmentBId = extractSegmentId(currentConfig.segmentB);
    const resolvedConfigResolution =
      resolutionSetting === "auto" ? "auto" : manualResolution;
    if (
      currentConfig.reference === reference &&
      currentConfig.track === track &&
      currentSegmentAId === segmentAId &&
      currentSegmentBId === segmentBId &&
      currentConfig.resolution === resolvedConfigResolution &&
      currentConfig.manualResolution === manualResolution &&
      currentConfig.autoResolution === autoResolution &&
      currentConfig.toggleColourScheme === toggleColourScheme
    ) {
      return;
    }
    props.handleViewUpdate(props.index, {
      ...props.viewConfig,
      config: {
        ...currentConfig,
        reference,
        track,
        segmentA: segmentAId,
        segmentB: segmentBId,
        resolution: resolvedConfigResolution,
        manualResolution,
        autoResolution,
        toggleColourScheme,
      },
    });
  }, [
    reference,
    track,
    segmentAId,
    segmentBId,
    resolutionSetting,
    manualResolution,
    autoResolution,
    toggleColourScheme,
  ]);

  const zoomRef = useRef<d3.ZoomBehavior<Element, unknown> | null>(null);
  const svgElRef = useRef<SVGSVGElement | null>(null);
  const automaticInitialRenderRef = useRef(false);

  type TileCoord = { x: number; y: number };

  interface RenderOptions {
    segmentA?: Segment;
    segmentB?: Segment;
    resolution?: number;
    locus?: {
      x: [number, number];
      y: [number, number];
    };
  }

  const determineAutoResolution = (xRange?: [number, number], yRange?: [number, number]) => {
    if (!xRange || !yRange) {
      return autoResolution || RESOLUTION_PRESETS[RESOLUTION_PRESETS.length - 1];
    }
    const xSpan = Math.max(1, xRange[1] - xRange[0]);
    const ySpan = Math.max(1, yRange[1] - yRange[0]);
    for (const candidate of RESOLUTION_PRESETS) {
      const tilesX = Math.max(1, Math.ceil(xSpan / (candidate * TILE_SIZE_BINS)));
      const tilesY = Math.max(1, Math.ceil(ySpan / (candidate * TILE_SIZE_BINS)));
      if (tilesX * tilesY <= AUTO_TILE_CEILING) {
        return candidate;
      }
    }
    return RESOLUTION_PRESETS[RESOLUTION_PRESETS.length - 1];
  };

  const resolveResolutionValue = (override?: number) => {
    if (typeof override === "number") return override;
    return resolutionSetting === "auto" ? autoResolution : manualResolution;
  };

  const ensureAutoResolution = (xRange: [number, number], yRange: [number, number]) => {
    const nextAuto = determineAutoResolution(xRange, yRange);
    setAutoResolution((prev) => (prev === nextAuto ? prev : nextAuto));
    return nextAuto;
  };

  const clearHeatmap = () => {
    const svg = d3.select(heatmapRef.current);
    svg.selectAll("*").remove();
  }

  const renderHeatmap = async (options: RenderOptions = {}) => {
    const targetSegmentA = options.segmentA ?? segmentA;
    const targetSegmentB = options.segmentB ?? segmentB;
    const targetResolution = resolveResolutionValue(options.resolution);
    const targetXLocus = options.locus?.x ?? xLocus;
    const targetYLocus = options.locus?.y ?? yLocus;

    if (!targetSegmentA || !targetSegmentB || !targetResolution) {
      console.warn("Missing data required to render heatmap tiles");
      return;
    }

    setLoading(true);
    try {
      const requiredTiles = getRequiredTiles(
        targetSegmentA,
        targetSegmentB,
        targetResolution,
        targetXLocus,
        targetYLocus,
      );
      const missingTiles = getMissingTiles(
        requiredTiles,
        targetSegmentA,
        targetSegmentB,
        targetResolution,
      );
      await fetchMissingTiles(
        missingTiles,
        targetSegmentA,
        targetSegmentB,
        targetResolution,
      );
      setRenderedSegmentA(targetSegmentA);
      setRenderedSegmentB(targetSegmentB);
      setRenderedResolution(targetResolution);
    } catch (error) {
      console.error("Failed to render heatmap tiles", error);
    } finally {
      setLoading(false);
    }
  };

  const handleZoomComplete = (nextXRange: [number, number], nextYRange: [number, number]) => {
    setXLocus(nextXRange);
    setYLocus(nextYRange);
    if (resolutionSetting !== "auto" || !segmentA || !segmentB) return;
    const nextResolution = ensureAutoResolution(nextXRange, nextYRange);
    renderHeatmap({
      resolution: nextResolution,
      locus: {
        x: nextXRange,
        y: nextYRange,
      },
    }).catch((error) => console.error("Failed to refresh tiles after zoom", error));
  };
  useEffect(() => {
    if (automaticInitialRenderRef.current) return;
    if (!reference || !track || !segmentA || !segmentB) return;
    automaticInitialRenderRef.current = true;
    const defaultXLocus: [number, number] = [0, segmentA.length];
    const defaultYLocus: [number, number] = [0, segmentB.length];
    setXLocus(defaultXLocus);
    setYLocus(defaultYLocus);
    const initialResolution =
      resolutionSetting === "auto"
        ? determineAutoResolution(defaultXLocus, defaultYLocus)
        : manualResolution;
    if (resolutionSetting === "auto") {
      setAutoResolution(initialResolution);
    }
    renderHeatmap({
      segmentA,
      segmentB,
      resolution: initialResolution,
      locus: { x: defaultXLocus, y: defaultYLocus },
    }).catch((error) => console.error("Failed to render heatmap on load", error));
  }, [reference, track, segmentA, segmentB, resolutionSetting, manualResolution]);

  const downloadPNG = () => {
    const canvas = canvasRef.current;
    const svgEl = svgElRef.current;
    if (!canvas || !svgEl) {
      console.warn("Canvas or SVG not ready");
      return;
    }
    const serializer = new XMLSerializer();
    let svgString = serializer.serializeToString(svgEl);
    if (!svgString.includes('xmlns="http://www.w3.org/2000/svg"')) {
      svgString = svgString.replace('<svg', '<svg xmlns="http://www.w3.org/2000/svg"');
    }
    const svgBlob = new Blob([svgString], { type: 'image/svg+xml;charset=utf-8' });
    const url = URL.createObjectURL(svgBlob);
    const img = new Image();
    img.onload = () => {
      const width = svgEl.clientWidth;
      const height = svgEl.clientHeight;
      const tmpCanvas = document.createElement('canvas');
      tmpCanvas.width = width;
      tmpCanvas.height = height;
      const ctx = tmpCanvas.getContext('2d');
      if (!ctx) return;
      
      ctx.fillStyle = '#ffffff';
      ctx.fillRect(0, 0, width, height);

      const style = getComputedStyle(canvas);
      const dx = parseInt(style.left, 10);
      const dy = parseInt(style.top, 10);
      ctx.drawImage(canvas, dx, dy);

      ctx.drawImage(img, 0, 0);
      URL.revokeObjectURL(url);
      const link = document.createElement('a');
      link.download = 'heatmap.png';
      link.href = tmpCanvas.toDataURL('image/png');
      link.click();
    };
    img.onerror = () => console.error('Failed to load SVG image for PNG export');
    img.src = url;
  };

  useEffect(() => {
    console.log(availableSegments);
  }, [availableSegments]);

    // Compute the required tiles based on current viewport
  const getRequiredTiles = (
    segA: Segment,
    segB: Segment,
    res: number,
    locusX: [number, number],
    locusY: [number, number],
  ): TileCoord[] => {
    if (!segA || !segB) return [];

    const binsX = Math.ceil(segA.length / res);
    const binsY = Math.ceil(segB.length / res);

    const tiles: TileCoord[] = [];
    const [xBpStart, xBpEnd] = locusX;
    const [yBpStart, yBpEnd] = locusY;

    const xBinStart = Math.floor(xBpStart / res);
    const xBinEnd = Math.ceil(xBpEnd / res);
    const yBinStart = Math.floor(yBpStart / res);
    const yBinEnd = Math.ceil(yBpEnd / res);

    const wx0 = Math.max(0, xBinStart);
    const wx1 = Math.min(binsX - 1, xBinEnd);
    const wy0 = Math.max(0, yBinStart);
    const wy1 = Math.min(binsY - 1, yBinEnd);

    if (wx0 > wx1 || wy0 > wy1) {
      return tiles;
    }

    const tileX0 = Math.floor(wx0 / TILE_SIZE_BINS);
    const tileX1 = Math.floor(wx1 / TILE_SIZE_BINS);
    const tileY0 = Math.floor(wy0 / TILE_SIZE_BINS);
    const tileY1 = Math.floor(wy1 / TILE_SIZE_BINS);

    for (let tx = tileX0; tx <= tileX1; tx++) {
      for (let ty = tileY0; ty <= tileY1; ty++) {
        tiles.push({ x: tx, y: ty });
      }
    }

    return tiles;
  };

  // Determine which tiles are missing from the cache
  const getMissingTiles = (
    requiredTiles: TileCoord[],
    segA: Segment,
    segB: Segment,
    res: number,
  ): TileCoord[] => {
    const missingTiles: TileCoord[] = [];
    if (!tileCache) return requiredTiles;
    for (const tile of requiredTiles) {
      const key = buildTileCacheKey({
        trackName: track,
        segmentAId: segA.id,
        segmentBId: segB.id,
        resolution: res,
        tileX: tile.x,
        tileY: tile.y,
      });
      if (!tileCache.has(key)) {
        missingTiles.push(tile);
      }
    }
    return missingTiles;
  };

  // Fetch missing tiles
  const fetchMissingTiles = async (
    missingTiles: TileCoord[],
    segA: Segment,
    segB: Segment,
    res: number,
  ) => {
    if (missingTiles.length === 0) return;

    const payload = {
      uuid: space.uuid,
      file_name: track,
      genome_path: reference,
      resolution: res,
      segment_1: segA.id,
      segment_2: segB.id,
      normalise,
      tiles: missingTiles,
    };

    try {
      const response = await fetch(buildApiUrl("/api/timge/heatmap_tiles/"), {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify(payload),
      });

      if (!response.ok) {
        throw new Error(`Tile request failed (${response.status})`);
      }

      const data = await response.json();
      if (data.status === "success" && data.tiles) {
        let cacheTouched = false;
        for (const tile of data.tiles as any[]) {
          const key = buildTileCacheKey({
            trackName: track,
            segmentAId: segA.id,
            segmentBId: segB.id,
            resolution: res,
            tileX: tile.x,
            tileY: tile.y,
          });
          tileCache?.set(key, tile.matrix ?? null);
          cacheTouched = true;
        }
        if (cacheTouched) {
          setTileCacheVersion((prev) => prev + 1);
        }
      } else {
        console.error("Failed to fetch tiles", data.message);
      }
    } catch (error) {
      console.error("Error fetching heatmap tiles", error);
      throw error;
    }
  };

  useEffect(() => {
    if (reference && track && props.dependencies) {
      const _segmentA = props.dependencies.segmentA as Segment | undefined;
      const _segmentB = props.dependencies.segmentB as Segment | undefined;
      if (_segmentA && _segmentB) {
        console.log("Dependencies changed", _segmentA, _segmentB);
        const nextXLocus: [number, number] = [0, _segmentA.length];
        const nextYLocus: [number, number] = [0, _segmentB.length];
        setSegmentAId(_segmentA.id);
        setSegmentBId(_segmentB.id);
        setXLocus(nextXLocus);
        setYLocus(nextYLocus);
        const nextResolution =
          resolutionSetting === "auto"
            ? ensureAutoResolution(nextXLocus, nextYLocus)
            : manualResolution;
        renderHeatmap({
          segmentA: _segmentA,
          segmentB: _segmentB,
          resolution: nextResolution,
          locus: { x: nextXLocus, y: nextYLocus },
        }).catch((error) => console.error("Failed to sync dependency tiles", error));
      }
    }
  }, [props.dependencies, reference, track, resolutionSetting, manualResolution]);

  return (
    <ParentView 
      viewConfig={props.viewConfig}
      index={props.index}
      userActions={{
        "Download PNG": downloadPNG,
        [props.viewConfig.config.isMinimised ? "Maximise" : "Minimise"]: () => {
          props.handleViewUpdate(props.index, {
            ...props.viewConfig,
            config: {
              ...props.viewConfig.config,
              isMinimised: !props.viewConfig.config.isMinimised,
            },
          });
        }
      }}
    >
        {!(reference && track) ? (
          // <Box sx={{ display: "flex", justifyContent: "center", alignItems: "center" }}>
          //   <Button variant="outlined" onClick={() => setShowTrackPicker(true)}>
          //     Select Tracks
          //   </Button>
          // </Box>
                      <Box
                      sx={{
                        display: "flex",
                        justifyContent: "center",
                        alignItems: "center",
                        width: "100%",
                        height: "100%",
                        flexDirection: "column",
                        flexGrow: 1,
                        padding: 2,
                      }}
                    >
                      <p>No Tracks Selected</p>
                      <Button
                        color="primary"
                        onClick={() => setShowTrackPicker(true)}
                        sx={{
                          marginTop: "10px",
                        }}
                      >
                        Select Tracks
                      </Button>
                    </Box>
        ) : (
        <Box className="flex flex-col gap-6 w-full"
            sx={{
                display: "flex",
                flexDirection: "column",
                justifyContent: "center",
                alignItems: "center",
                width: "100%",
                height: "100%",      
                flexGrow: 1,  
            }}
        >
            <Card
            variant="outlined"
            className="w-full"
            sx={{
                borderRadius: "none",
                boxShadow: "none",
                width: "100%",
                backgroundColor: "#f3f3f3",
                border: "1px solid #bfbfbf",
            }}
            >
              <Box
                sx={{
                  display: "flex",
                  justifyContent: "space-between",
                  alignItems: "center",
                  padding: "10px",
                }}
              >
                <Box className="flex flex-wrap gap-4 items-center"
                  sx={{
                    display: "flex",
                    flexWrap: "wrap",
                  }}
                >
                <Box
                    sx={{
                      display: "flex",
                      flexDirection: "row",
                      alignItems: "center",
                      justifyContent: "left",
                      gap: "10px",
                    }}
                  >
                    <Typography 
                      sx={{
                        fontSize: "0.8em",
                      }}
                    >
                    Segment A:
                    </Typography>
                    <Select
                    onChange={(e, value) => {
                        dispatch(updateViewConfig({
                            uuid: props.viewConfig.uuid,
                            config: {
                                ...props.viewConfig.config,
                                segmentA: value,
                            }
                        }))
                        setSegmentAId(value ?? null);
                    }}
                    value={segmentAId ?? null}
                    placeholder="Select segment A"
                    sx={{
                      boxShadow: "none",
                      fontSize: "0.8em",
                    }}
                    >
                    {Array.from(availableSegments).map((segment) => (
                        <Option key={segment.id} value={segment.id}
                          sx={{
                            fontSize: "0.8em",
                          }}
                        >
                        {segment.id}
                        </Option>
                    ))}
                    </Select>
                  </Box>
                  <Box
                    sx={{
                      display: "flex",
                      flexDirection: "row",
                      alignItems: "center",
                      justifyContent: "left",
                      gap: "10px",
                    }}
                  >
                    <Typography 
                      sx={{
                        fontSize: "0.8em",
                      }}
                    >
                    Segment B:
                    </Typography>
                    <Select
                    onChange={(e, value) => {
                        dispatch(updateViewConfig({
                            uuid: props.viewConfig.uuid,
                            config: {
                                ...props.viewConfig.config,
                                segmentB: value,
                            }
                        }))
                        setSegmentBId(value ?? null);
                    }}
                    value={segmentBId ?? null}
                    placeholder="Select segment B"
                    sx={{
                      boxShadow: "none",
                      fontSize: "0.8em",
                    }}
                    >
                    {Array.from(availableSegments).map((segment) => (
                        <Option key={segment.id} value={segment.id}
                          sx={{
                            fontSize: "0.8em",
                          }}
                          >
                        {segment.id}
                        </Option>
                    ))}
                    </Select>
                    </Box>
                    <Box
                    sx={{
                      display: "flex",
                      flexDirection: "row",
                      alignItems: "center",
                      justifyContent: "left",
                      gap: "10px",
                    }}
                  >
                    <Typography 
                      sx={{
                        fontSize: "0.8em",
                      }}
                    >
                    Resolution (bp):
                    </Typography>
                    <Select
                      value={resolutionSetting === "auto" ? "auto" : manualResolution}
                      placeholder="Select resolution"
                      onChange={(e, value) => {
                        if (!value) return;
                        if (value === "auto") {
                          setResolutionSetting("auto");
                          if (segmentA && segmentB) {
                            const nextResolution = ensureAutoResolution(xLocus, yLocus);
                            renderHeatmap({
                              resolution: nextResolution,
                              locus: { x: xLocus, y: yLocus },
                            }).catch((error) =>
                              console.error("Failed to render heatmap in auto mode", error),
                            );
                          }
                          return;
                        }
                        setResolutionSetting(value);
                        setManualResolution(value);
                      }}
                      sx={{
                        boxShadow: "none",
                        fontSize: "0.8em",
                      }}
                    >
                      <Option value="auto" sx={{ fontSize: "0.8em" }}>
                        Auto
                      </Option>
                      {RESOLUTION_PRESETS.map((resOption) => (
                        <Option key={resOption} value={resOption} sx={{ fontSize: "0.8em" }}>
                          {resOption}
                        </Option>
                      ))}
                    </Select>
                  </Box>
                  <Box
                    sx={{
                      display: "flex",
                      flexDirection: "row",
                      alignItems: "center",
                      justifyContent: "left",
                      gap: "10px",
                    }}
                  >
                    <Typography
                      sx={{
                        fontSize: "0.8em",
                      }}
                    >
                    Toggle colour scheme:
                    </Typography>
                    <Checkbox
                    checked={toggleColourScheme}
                    onChange={(e) => {
                        setToggleColourScheme(e.target.checked);
                    }}
                    />
                  </Box>
                  <Box
                    sx={{
                      display: "flex",
                      flexDirection: "row",
                      alignItems: "center",
                      justifyContent: "left",
                      gap: "10px",
                    }}
                  >
                    <Typography 
                      sx={{
                        fontSize: "0.8em",
                      }}
                    >
                    Normalise:
                    </Typography>
                    <Checkbox
                    checked={normalise}
                    onChange={(e) => {
                        setNormalise(e.target.checked);
                    }}
                    />
                  </Box>
                  <Box
                    sx={{
                      display: "flex",
                      flexDirection: "row",
                      alignItems: "center",
                      justifyContent: "left",
                      gap: "10px",
                    }}
                  >
                      <Typography 
                        sx={{
                          fontSize: "0.8em",
                        }}
                      >
                      Negative Strand:
                      </Typography>
                      <Checkbox
                      checked={negativeStrand}
                      onChange={(e) => {
                          setNegativeStrand(e.target.checked);
                          dispatch(updateViewConfig({
                            uuid: props.viewConfig.uuid,
                            config: {
                                ...props.viewConfig.config,
                                negativeStrand: e.target.checked,
                            }
                        }))
                      }}
                      />
                    </Box>
                    <Box
                    sx={{
                      display: "flex",
                      flexDirection: "row",
                      alignItems: "center",
                      justifyContent: "left",
                      gap: "10px",
                    }}
                  >
                    <Typography
                      sx={{
                        fontSize: "0.8em",
                      }}
                    >Locus:</Typography>
                    <Box display="flex" gap={2}>
                      {/* <input
                        type="text"
                        value={xLocus.join(" - ")}
                        onChange={(e) => {
                          const [start, end] = e.target.value.split("-").map(s => +s.trim());
                          setXLocus([start, end]);
                        }}
                        onBlur={() => {
                          if (zoomToLocusRef.current) {
                            zoomToLocusRef.current(xLocus[0], xLocus[1], yLocus[0], yLocus[1]);
                          }
                        }}
                      />
                      <input
                        type="text"
                        value={yLocus.join(" - ")}
                        onChange={(e) => {
                          const [start, end] = e.target.value.split("-").map(s => +s.trim());
                          setYLocus([start, end]);
                        }}
                        onBlur={() => {
                          if (zoomToLocusRef.current) {
                            zoomToLocusRef.current(xLocus[0], xLocus[1], yLocus[0], yLocus[1]);
                          }
                        }}
                      /> */}
                      <Typography
                        sx={{
                          fontSize: "0.8em",
                        }}
                      >{xLocus[0]} - {xLocus[1]},</Typography>
                      <Typography
                        sx={{
                          fontSize: "0.8em",
                        }}
                      >{yLocus[0]} - {yLocus[1]}</Typography>
                    </Box>
                    </Box>
                    <Button variant="solid" color="primary" 
                      onClick={() => {
                        const nextXLocus: [number, number] = [0, segmentA?.length ?? 0];
                        const nextYLocus: [number, number] = [0, segmentB?.length ?? 0];
                        setXLocus(nextXLocus);
                        setYLocus(nextYLocus);
                        clearHeatmap();
                        const nextResolution =
                          resolutionSetting === "auto"
                            ? ensureAutoResolution(nextXLocus, nextYLocus)
                            : manualResolution;
                        renderHeatmap({
                          resolution: nextResolution,
                          locus: { x: nextXLocus, y: nextYLocus },
                        });
                      }}
                      sx={{
                        fontSize: "0.8em",
                      }}
                    >
                    Render
                    </Button>
                    <Button
                      variant="outlined"
                      color="neutral"
                      onClick={() => {
                        if (zoomRef.current && svgElRef.current) {
                          d3.select(svgElRef.current)
                            .transition()
                            .duration(500)
                            .call(zoomRef.current.transform, d3.zoomIdentity);
                        } else {
                          console.warn("Zoom or SVG not initialized");
                        }
                      }}
                      sx={{
                        fontSize: "0.8em",
                      }}
                    >
                      Reset Zoom
                    </Button>
                </Box>
                </Box>
            </Card>
            <Box
              sx={{
                position: "relative",
                height: "100%",
                width: "100%",
                display: "flex",
                justifyContent: "center",
                alignItems: "center",
                marginBottom: "20px",
                padding: "10px",
              }}
            >
              {loading && (
                <Box
                  sx={{
                    position: "absolute",
                    inset: 0,
                    display: "flex",
                    alignItems: "center",
                    justifyContent: "center",
                    backgroundColor: "rgba(255, 255, 255, 0.75)",
                    zIndex: 2,
                    borderRadius: "8px",
                  }}
                >
                  <CircularProgress size="lg" variant="solid"/>
                </Box>
              )}
              <Box
                sx={{
                  display: "flex",
                  width: "100%",
                  height: "100%",
                  justifyContent: "center",
                  alignItems: "center",
                  opacity: loading ? 0.4 : 1,
                  transition: "opacity 0.2s ease",
                }}
              >
                <CanvasHeatmap
                  setCanvasRef={(el) => {
                    canvasRef.current = el;
                  }}
                  setZoomRef={(zoom, svgEl) => {
                    zoomRef.current = zoom;
                    svgElRef.current = svgEl;
                  }}
                  zoomToLocusRef={zoomToLocusRef}
                  onLocusChange={(xRange, yRange) => {
                    setXLocus(xRange);
                    setYLocus(yRange);
                  }}
                  onZoomEnd={handleZoomComplete}
                  title={`${track.split("/").pop()} (${renderedResolution}nt)`}
                  trackName={track}
                  segmentA={renderedSegmentA}
                  segmentB={renderedSegmentB}
                  resolution={renderedResolution}
                  tileCacheVersion={tileCacheVersion}
                  toggleColourScheme={toggleColourScheme}
                  isMinimised={props.viewConfig.config.isMinimised}
                />
              </Box>
            </Box>
        </Box>
            )
        }
        {showTrackPicker && (
          <TrackSelector
            onClose={() => setShowTrackPicker(false)}
            onConfirm={handleTrackConfirm}
          />
        )}
    </ParentView>

  );
};

export default MapView;
