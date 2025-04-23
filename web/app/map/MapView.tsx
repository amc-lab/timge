"use client";
import React, { useEffect, useRef, useState } from "react";
import * as d3 from "d3";
import ParentView from "@/components/ParentView";
import { Box, Button, Card, Checkbox, CircularProgress, Dropdown, LinearProgress, Option, Select, Typography } from "@mui/joy";
import { useAppDispatch, useAppSelector } from "@/store/hooks";
import { View } from "@/store/features/views/types";
import TrackSelector from "./components/TrackSelector";

interface MapViewProps {
  trackFiles: any[];
  viewConfig: View;
  handleViewUpdate: (index, viewState: View) => void;
  index: number;
  crossViewActionHandler?: any;
  dependencies: any;
  addConnection?: any;
  removeConnection?: any;
  // createdViews: Set<any>;
}

const MapView = (props: MapViewProps) => {
  const [reference, setReference] = useState(props.viewConfig.config.reference);
  const [track, setTrack] = useState(props.viewConfig.config.track);
  const [segmentA, setSegmentA] = useState(props.viewConfig.config.segmentA);
  const [segmentB, setSegmentB] = useState(props.viewConfig.config.segmentB);
  const [resolution, setResolution] = useState(props.viewConfig.config.resolution);
  const [availableSegments, setAvailableSegments] = useState([]);
  const [showGridlines, setShowGridlines] = useState(false);
  const [toggleColourScheme, setToggleColourScheme] = useState(true); // default to white-red
  const [normalise, setNormalise] = useState(false);
  const [loading, setLoading] = useState(false);

  const heatmapRef = useRef(null);
  const space = useAppSelector((state) => state.space);

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
      },
    });
  }
  , [reference, track, segmentA, segmentB, resolution]);

  const zoomRef = useRef<d3.ZoomBehavior<Element, unknown> | null>(null);

  const drawHeatmap = (matrix: number[][], segmentA: string, segmentB: string) => {
    setLoading(true);
    const svg = d3.select(heatmapRef.current);
    svg.selectAll("*").remove();

    const maxWidth = window.innerWidth * 0.8 * (props.viewConfig.config.isMinimised ? 0.5 : 1) - 150;
    const maxHeight = window.innerHeight * 0.8;
  
    const numRows = matrix.length;
    const numCols = matrix[0].length;
  
    const cellSize = Math.min(maxWidth / numRows, maxHeight / numCols);
  
    const margin = { top: 20, right: 80, bottom: 70, left: 70 };
    const width = numRows * cellSize;
    const height = numCols * cellSize;
  
    svg
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
      .attr("viewBox", [0, 0, width + margin.left + margin.right, height + margin.top + margin.bottom]);
  
    const g = svg.append("g")
      .attr("transform", `translate(${margin.left},${margin.top})`);
  
    const flatValues = matrix.flat();
    const colorScale = d3.scaleSequential()
      .domain([0, d3.max(flatValues)!])
      .interpolator(toggleColourScheme ? d3.interpolateOrRd : d3.interpolateGnBu);
  
    const xScale = d3.scaleLinear()
      .domain([0, numRows])
      .range([0, width]);
  
    const yScale = d3.scaleLinear()
      .domain([0, numCols])
      .range([height, 0]);
  
    const xAxis = d3.axisBottom(xScale)
      .ticks(Math.min(numRows, 10))
      .tickFormat(d => `${d * resolution}`);
  
    const yAxis = d3.axisLeft(yScale)
      .ticks(Math.min(numCols, 10))
      .tickFormat(d => `${d * resolution}`);
  
    const gx = svg.append("g")
      .attr("transform", `translate(${margin.left}, ${height + margin.top})`)
      .attr("class", "x-axis")
      .call(xAxis);
  
    const gy = svg.append("g")
      .attr("transform", `translate(${margin.left}, ${margin.top})`)
      .attr("class", "y-axis")
      .call(yAxis);
  
    const heatmap = g.append("g").attr("class", "heatmap");
  
    heatmap.selectAll("g")
      .data(matrix)
      .join("g")
      .attr("transform", (_, i) => `translate(${(numRows - 1 - i) * cellSize}, 0)`)
      .selectAll("rect")
      .data((d, i) => d.map((val, j) => ({ val, j })))
      .join("rect")
      .attr("y", d => (d.j) * cellSize)
      .attr("width", cellSize)
      .attr("height", cellSize)
      .attr("fill", d => colorScale(d.val))
      .attr("stroke", showGridlines ? "#ccc" : "none");
  
    svg.append("text")
      .attr("x", margin.left + width / 2)
      .attr("y", height + margin.top + 35)
      .attr("text-anchor", "middle")
      .attr("font-size", "12px")
      .text(segmentA);
  
    svg.append("text")
      .attr("transform", `rotate(-90)`)
      .attr("x", -margin.top - height / 2)
      .attr("y", 12)
      .attr("text-anchor", "middle")
      .attr("font-size", "12px")
      .text(segmentB);
  
    const legendHeight = 200;
    const legendWidth = 15;
  
    const legendScale = d3.scaleLinear()
      .domain(colorScale.domain())
      .range([legendHeight, 0]);
  
    const legendAxis = d3.axisRight(legendScale)
      .ticks(5);
  
    const legend = svg.append("g")
      .attr("class", "legend")
      .attr("transform", `translate(${margin.left + width + 20}, ${margin.top})`);
  
    const legendGradientId = "legend-gradient";
    setLoading(false);
  
    svg.select("defs").remove();
    const defs = svg.append("defs");
    const linearGradient = defs.append("linearGradient")
      .attr("id", legendGradientId)
      .attr("x1", "0%")
      .attr("y1", "100%")
      .attr("x2", "0%")
      .attr("y2", "0%");
  
    const numStops = 10;
    const step = 1 / (numStops - 1);
  
    d3.range(numStops).forEach(i => {
      linearGradient.append("stop")
        .attr("offset", `${i * step * 100}%`)
        .attr("stop-color", colorScale(colorScale.domain()[0] + i * step * (colorScale.domain()[1] - colorScale.domain()[0])));
    });
  
    legend.append("rect")
      .attr("width", legendWidth)
      .attr("height", legendHeight)
      .style("fill", `url(#${legendGradientId})`);
  
    legend.append("g")
      .attr("transform", `translate(${legendWidth}, 0)`)
      .call(legendAxis);
    
    svg.style("shape-rendering", "crispEdges");

    zoomRef.current = d3.zoom()
      .scaleExtent([1, 10])
      .translateExtent([[-100, -100], [width + 100, height + 100]])
      .on("zoom", zoomed);
  
    svg.call(zoomRef.current);
  
    function zoomed({ transform }) {
      heatmap.attr("transform", transform);
      const zx = transform.rescaleX(xScale);
      const zy = transform.rescaleY(yScale);
      gx.call(xAxis.scale(zx));
      gy.call(yAxis.scale(zy));
    }

  };

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
          drawHeatmap(data.matrix, _segmentA ? _segmentA : segmentA, _segmentB ? _segmentB : segmentB);
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
        setSegmentA(props.dependencies.segmentA);
        setSegmentB(props.dependencies.segmentB);
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
          <Box sx={{ display: "flex", justifyContent: "center", alignItems: "center" }}>
            <Button variant="outlined" onClick={() => setShowTrackPicker(true)}>
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
                        setSegmentA(value)
                        props.handleViewUpdate(props.index, {
                            ...props.viewConfig,
                            config: {
                                ...props.viewConfig.config,
                                segmentA: value,
                            },
                        })
                    }}
                    value={segmentA}
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
                        setSegmentB(value)
                        props.handleViewUpdate(props.index, {
                            ...props.viewConfig,
                            config: {
                                ...props.viewConfig.config,
                                segmentB: value,
                            },
                        })
                    }}
                    value={segmentB}
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
                    {/* <Option value={1}>1</Option> */}
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
              {loading && (
                <Box
                  sx={{
                    display: "flex",
                    justifyContent: "center",
                    alignItems: "center",
                  }}
                >
                  <CircularProgress
                  />
                </Box>
            )} {
              <svg ref={heatmapRef}></svg>
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
