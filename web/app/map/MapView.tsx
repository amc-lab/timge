"use client";
import React, { useEffect, useRef, useState } from "react";
import * as d3 from "d3";
import ParentView from "@/components/ParentView";
import { Box, Button, Card, Checkbox, CircularProgress, Dropdown, LinearProgress, Option, Select, Typography } from "@mui/joy";
import { useAppDispatch, useAppSelector } from "@/store/hooks";
import { View } from "@/store/features/views/types";
import TrackSelector from "./components/TrackSelector";
import CanvasHeatmap from "./components/CanvasHeatmap";
import { getViewConfig, updateViewConfig } from "@/store/features/space/spaceSlice";

interface MapViewProps {
  trackFiles: any[];
  viewConfig: View;
  handleViewUpdate: (index, viewState: View) => void;
  index: number;
  crossViewActionHandler?: any;
  dependencies: any;
  addConnection?: any;
  removeConnection?: any;
}

const MapView = (props: MapViewProps) => {
  const [reference, setReference] = useState(props.viewConfig.config.reference);
  const [track, setTrack] = useState(props.viewConfig.config.track);
  // const [segmentA, setSegmentA] = useState(props.viewConfig.config.segmentA);
  // const [segmentB, setSegmentB] = useState(props.viewConfig.config.segmentB);
  const [resolution, setResolution] = useState(props.viewConfig.config.resolution);
  const [availableSegments, setAvailableSegments] = useState([]);
  const [showGridlines, setShowGridlines] = useState(false);
  const [toggleColourScheme, setToggleColourScheme] = useState(props.viewConfig.config.toggleColourScheme || false);
  const [normalise, setNormalise] = useState(false);
  const [negativeStrand, setNegativeStrand] = useState(props.viewConfig.config.negativeStrand || false);
  const [loading, setLoading] = useState(false);
  const [matrix, setMatrix] = useState<number[][]>([]);
  const canvasRef = useRef<HTMLCanvasElement | null>(null);

  const heatmapRef = useRef(null);
  const space = useAppSelector((state) => state.space);
  const dispatch = useAppDispatch();

  // const segmentA = space.views[props.index].config.segmentA || props.viewConfig.config.segmentA;
  // const segmentB = space.views[props.index].config.segmentB || props.viewConfig.config.segmentB;

  const [segmentA, setSegmentA] = useState(space.views[props.index].config.segmentA || props.viewConfig.config.segmentA);
  const [segmentB, setSegmentB] = useState(space.views[props.index].config.segmentB || props.viewConfig.config.segmentB);


  const [renderedSegmentA, setRenderedSegmentA] = useState(segmentA);
  const [renderedSegmentB, setRenderedSegmentB] = useState(segmentB);
  const [renderedResolution, setRenderedResolution] = useState(resolution);

  const [showTrackPicker, setShowTrackPicker] = useState(false);

  const [xLocus, setXLocus] = useState<[number, number]>([0, 0]);
  const [yLocus, setYLocus] = useState<[number, number]>([0, 0]);
  const zoomToLocusRef = useRef<((x0: number, x1: number, y0: number, y1: number) => void) | null>(null);


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
      const host = process.env.NEXT_PUBLIC_DJANGO_HOST;
      const queryParams = new URLSearchParams({
        uuid: space.uuid,
        file_name: track,
        genome_path: reference,
      }).toString();
        const response = await fetch(`${host}/api/timge/get_segments/?${queryParams}`, {
            method: "GET",
            headers: {
            "Content-Type": "application/json",
            },
        });
      const data = await response.json();
      if (data.status === "success" && data.segments) {
        setAvailableSegments(data.segments);
      } else {
        console.error("Failed to fetch segments", data.message);
      }
    };

    if (reference && track) {
      fetchData();
    }
  }, [reference, track]);

  useEffect(() => {
    props.handleViewUpdate(props.index, {
      ...props.viewConfig,
      config: {
        ...props.viewConfig.config,
        reference: reference,
        track: track,
        segmentA: segmentA,
        segmentB: segmentB,
        resolution: resolution,
        toggleColourScheme: toggleColourScheme,
      },
    });
  }
  , [reference, track, segmentA, segmentB, resolution, toggleColourScheme]);

  const zoomRef = useRef<d3.ZoomBehavior<Element, unknown> | null>(null);
  const svgElRef = useRef<SVGSVGElement | null>(null);


  const clearHeatmap = () => {
    const svg = d3.select(heatmapRef.current);
    svg.selectAll("*").remove();
  }
  
  const renderHeatmap = (_segmentA?, _segmentB?) => {
    setLoading(true);
    console.log("Rendering heatmap", _segmentA, _segmentB);
    const host = process.env.NEXT_PUBLIC_DJANGO_HOST;
    fetch(`${host}/api/timge/heatmap/`, {
      method: "POST",
      headers: {
        "Content-Type": "application/json"
      },
      body: JSON.stringify({
        uuid: space.uuid,
        file_name: track,
        genome_path: reference,
        resolution: resolution,
        segment_1: _segmentA ? _segmentA : segmentA,
        segment_2: _segmentB ? _segmentB : segmentB,
        normalise: normalise,
      })
    })
      .then(response => response.json())
      .then(data => {
        if (data.status === "success") {
          if (negativeStrand) {
            const numRows = data.matrix.length;
            const numCols = data.matrix[0]?.length || 0;
            setMatrix(
              data.matrix.map((row, rowIndex) =>
                row.map((_, colIndex) =>
                  data.matrix[numRows - 1 - rowIndex][numCols - 1 - colIndex]
                )
              )
            );
          } else {
            setMatrix(data.matrix);
          }
          setRenderedSegmentA(_segmentA ? _segmentA : segmentA); 
          setRenderedSegmentB(_segmentB ? _segmentB : segmentB);
          setRenderedResolution(resolution);
        } else {
          console.error("Failed to generate heatmap", data.message);
        }
      }
    )
  }

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
    if (reference && track) {
      if (props.dependencies) {
        const _segmentA = props.dependencies.segmentA;
        const _segmentB = props.dependencies.segmentB;
        console.log("Dependencies changed", _segmentA, _segmentB);
        setSegmentA(_segmentA);
        setSegmentB(_segmentB);
        renderHeatmap(_segmentA, _segmentB);
      }
    }
  }
  , [props.dependencies]);

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
                        setSegmentA(value);
                    }}
                    value={segmentA}
                    placeholder="Select segment A"
                    sx={{
                      boxShadow: "none",
                      fontSize: "0.8em",
                    }}
                    >
                    {availableSegments.map((segment, index) => (
                        <Option key={index} value={segment}
                          sx={{
                            fontSize: "0.8em",
                          }}
                        >
                        {segment}
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
                        setSegmentB(value);
                    }}
                    value={segmentB}
                    placeholder="Select segment B"
                    sx={{
                      boxShadow: "none",
                      fontSize: "0.8em",
                    }}
                    >
                    {availableSegments.map((segment, index) => (
                        <Option key={index} value={segment}
                          sx={{
                            fontSize: "0.8em",
                          }}
                          >
                        {segment}
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
                    defaultValue={resolution}
                    placeholder="Select resolution"
                    onChange={(e, value) => {
                        setResolution(value)
                        props.handleViewUpdate(props.index, {
                            ...props.viewConfig,
                            config: {
                                ...props.viewConfig.config,
                                resolution: value,
                            },
                        })
                    }}
                    sx={{
                      boxShadow: "none",
                      fontSize: "0.8em",
                    }}
                    >
                    <Option value={1} sx={{fontSize: "0.8em"}}>1</Option>
                    <Option value={5} sx={{fontSize: "0.8em"}}>5</Option>
                    <Option value={10} sx={{fontSize: "0.8em"}}>10</Option>
                    <Option value={25} sx={{fontSize: "0.8em"}}>25</Option>
                    <Option value={50} sx={{fontSize: "0.8em"}}>50</Option>
                    <Option value={100} sx={{fontSize: "0.8em"}}>100</Option>
                    <Option value={200} sx={{fontSize: "0.8em"}}>200</Option>
                    <Option value={500} sx={{fontSize: "0.8em"}}>500</Option>
                    <Option value={1000} sx={{fontSize: "0.8em"}}>1000</Option>
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
                    Show gridlines:
                    </Typography>
                    <Checkbox
                    checked={showGridlines}
                    onChange={(e) => {
                        setShowGridlines(e.target.checked);
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
                        clearHeatmap();
                        renderHeatmap(null, null)
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
                height: '100%',
                width: '100%',
                display: 'flex',
                justifyContent: 'center',
                alignItems: 'center',
                marginBottom: '20px',
                padding: '10px',
              }}
            >
               {
              <Box
                sx={{
                  display: "flex",
                  width: "100%",
                  height: "100%",
                  justifyContent: "center",
                  alignItems: "center",
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
                title={`${track.split("/").pop()} (${renderedResolution}nt)`}
                matrix={matrix}
                segmentA={renderedSegmentA}
                segmentB={renderedSegmentB}
                resolution={renderedResolution}
                toggleColourScheme={toggleColourScheme}
                showGridlines={showGridlines}
                isMinimised={props.viewConfig.config.isMinimised}
              />
              </Box>
            }
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
