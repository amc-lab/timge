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
  console.log(props.viewConfig);
  const [reference, setReference] = useState(props.viewConfig.config.reference);
  const [track, setTrack] = useState(props.viewConfig.config.track);
  // const [segmentA, setSegmentA] = useState(props.viewConfig.config.segmentA);
  // const [segmentB, setSegmentB] = useState(props.viewConfig.config.segmentB);
  const [resolution, setResolution] = useState(props.viewConfig.config.resolution);
  const [availableSegments, setAvailableSegments] = useState([]);
  const [showGridlines, setShowGridlines] = useState(false);
  const [toggleColourScheme, setToggleColourScheme] = useState(props.viewConfig.config.toggleColourScheme || false);
  const [normalise, setNormalise] = useState(false);
  const [loading, setLoading] = useState(false);
  const [matrix, setMatrix] = useState<number[][]>([]);

  const heatmapRef = useRef(null);
  const space = useAppSelector((state) => state.space);
  const dispatch = useAppDispatch();

  const segmentA = space.views[props.index].config.segmentA || props.viewConfig.config.segmentA;
  const segmentB = space.views[props.index].config.segmentB || props.viewConfig.config.segmentB;

  const [renderedSegmentA, setRenderedSegmentA] = useState(segmentA);
  const [renderedSegmentB, setRenderedSegmentB] = useState(segmentB);
  const [renderedResolution, setRenderedResolution] = useState(resolution);

  const [showTrackPicker, setShowTrackPicker] = useState(false);

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

  // useEffect(() => {
  //   if (props.viewConfig.config.segmentA && props.viewConfig.config.segmentB && props.viewConfig.config.resolution) {
  //     renderHeatmap();
  //   }
  // }, []);

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
          setMatrix(data.matrix);
          setRenderedSegmentA(_segmentA ? _segmentA : segmentA); 
          setRenderedSegmentB(_segmentB ? _segmentB : segmentB);
          setRenderedResolution(resolution);
        } else {
          console.error("Failed to generate heatmap", data.message);
        }
      }
    )
  }

  useEffect(() => {
    if (reference && track) {
      if (props.dependencies) {
        const _segmentA = props.dependencies.segmentA;
        const _segmentB = props.dependencies.segmentB;
        console.log("Dependencies changed", _segmentA, _segmentB);
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
        "Download PNG": () => {
              const d3ToPng = require('d3-svg-to-png');
              d3ToPng(heatmapRef.current, 'output', {
                scale: 1,
                format: 'png',
                quality: 1,
                download: true,
                ignore: '.ignored',
                background: 'white'
              });
            },
            [props.viewConfig.config.isMinimised ? "Maximise" : "Minimise"]: () => {
              clearHeatmap();
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
            }}
        >
            <Card
            variant="outlined"
            className="w-full"
            sx={{
                borderRadius: "none",
                boxShadow: "none",
                width: "100%",
            }}
            >
                <Box className="flex flex-wrap gap-4 items-center">
                    <Typography fontSize="md">
                    Segment A:
                    </Typography>
                    <Select
                    onChange={(e, value) => {
                        // setSegmentA(value)
                        dispatch(updateViewConfig({
                            uuid: props.viewConfig.uuid,
                            config: {
                                ...props.viewConfig.config,
                                segmentA: value,
                            }
                        }))
                        // props.handleViewUpdate(props.index, {
                        //     ...props.viewConfig,
                        //     config: {
                        //         ...props.viewConfig.config,
                        //         segmentA: value,
                        //     },
                        // })
                    }}
                    value={props.viewConfig.config.segmentA}
                    placeholder="Select segment A"
                    sx={{
                      boxShadow: "none",
                    }}
                    >
                    {availableSegments.map((segment, index) => (
                        <Option key={index} value={segment}>
                        {segment}
                        </Option>
                    ))}
                    </Select>

                    <Typography fontSize="md">
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
                    }}
                    value={props.viewConfig.config.segmentB}
                    placeholder="Select segment B"
                    sx={{
                      boxShadow: "none",
                    }}
                    >
                    {availableSegments.map((segment, index) => (
                        <Option key={index} value={segment}>
                        {segment}
                        </Option>
                    ))}
                    </Select>

                    <Typography fontSize="md">
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
                    }}
                    >
                    <Option value={1}>1</Option>
                    <Option value={5}>5</Option>
                    <Option value={10}>10</Option>
                    <Option value={25}>25</Option>
                    <Option value={50}>50</Option>
                    <Option value={100}>100</Option>
                    <Option value={200}>200</Option>
                    <Option value={500}>500</Option>
                    <Option value={1000}>1000</Option>
                    </Select>
                    <Typography fontSize="md">
                    Show gridlines:
                    </Typography>
                    <Checkbox
                    checked={showGridlines}
                    onChange={(e) => {
                        setShowGridlines(e.target.checked);
                    }}
                    />
                    <Typography fontSize="md">
                    Toggle colour scheme:
                    </Typography>
                    <Checkbox
                    checked={toggleColourScheme}
                    onChange={(e) => {
                        setToggleColourScheme(e.target.checked);
                    }}
                    />
                    <Typography fontSize="md">
                    Normalise:
                    </Typography>
                    <Checkbox
                    checked={normalise}
                    onChange={(e) => {
                        setNormalise(e.target.checked);
                    }}
                    />
                    <Button variant="solid" color="primary" onClick={() => {
                      clearHeatmap();
                      renderHeatmap(null, null)
                    }
                    }>
                    Render
                    </Button>
                    <Button
                      variant="outlined"
                      color="neutral"
                      onClick={() => {
                        if (zoomRef.current) {
                          const svg = d3.select(heatmapRef.current);
                          svg.transition()
                            .duration(500)
                            .call(zoomRef.current.transform, d3.zoomIdentity);
                        } else {
                          console.warn("Zoom behavior not initialized yet");
                        }
                      }}
                    >
                      Reset Zoom
                    </Button>
                </Box>
            </Card>
            <Box
              sx={{
                maxWidth: '90vw',
                maxHeight: '90vh',
                overflow: 'auto',
                // border: '1px solid #ccc',
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
                }}
                >
              <CanvasHeatmap
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
