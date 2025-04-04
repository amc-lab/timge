"use client";
import React, { useEffect, useRef, useState } from "react";
import * as d3 from "d3";
import { View } from "../types/state";
import ParentView from "@/components/ParentView";
import { Box, Button, Card, Dropdown, Option, Select, Typography } from "@mui/joy";

interface MapViewProps {
  trackFiles: any[];
  viewConfig: View;
  handleViewUpdate: (index, viewState: View) => void;
  index: number;
}

const MapView = (props: MapViewProps) => {
  const [reference, setReference] = useState("");
  const [track, setTrack] = useState("");
  const [segmentA, setSegmentA] = useState("");
  const [segmentB, setSegmentB] = useState("");
  const [resolution, setResolution] = useState(25);
  const [availableSegments, setAvailableSegments] = useState([]);

  const heatmapRef = useRef(null);

  useEffect(() => {
    const fetchData = async () => {
      const host = process.env.NEXT_PUBLIC_DJANGO_HOST;
      const queryParams = new URLSearchParams({
        uuid: props.viewConfig.uuid,
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

  let zoom;

  const drawHeatmap = (matrix: number[][]) => {
    const svg = d3.select(heatmapRef.current);
    svg.selectAll("*").remove(); // Clear previous drawing
  
    const cellSize = 20;
    const numRows = matrix.length;
    const numCols = matrix[0].length;
  
    const margin = { top: 20, right: 20, bottom: 40, left: 40 };
    const width = numCols * cellSize;
    const height = numRows * cellSize;
  
    svg
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
      .attr("viewBox", [0, 0, width + margin.left + margin.right, height + margin.top + margin.bottom]);
  
    const g = svg.append("g")
      .attr("transform", `translate(${margin.left},${margin.top})`);
  
    const flatValues = matrix.flat();
    const colorScale = d3.scaleSequential()
      .domain([0, d3.max(flatValues)!])
      .interpolator(d3.interpolateOrRd);
  
    const xScale = d3.scaleLinear()
      .domain([0, numCols])
      .range([0, width]);
  
    const yScale = d3.scaleLinear()
      .domain([0, numRows])
      .range([0, height]);
  
    const xAxis = d3.axisBottom(xScale)
      .ticks(Math.min(numCols, 10))
      .tickFormat(d => `${d * resolution}`);
  
    const yAxis = d3.axisLeft(yScale)
      .ticks(Math.min(numRows, 10))
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
  
    // Draw heatmap cells
    heatmap.selectAll("g")
      .data(matrix)
      .join("g")
      .attr("transform", (_, i) => `translate(0, ${i * cellSize})`)
      .selectAll("rect")
      .data(d => d)
      .join("rect")
      .attr("x", (_, j) => j * cellSize)
      .attr("width", cellSize)
      .attr("height", cellSize)
      .attr("fill", d => colorScale(d))
      .attr("stroke", "#ccc");
  
    // Axis labels
    svg.append("text")
      .attr("x", margin.left + width / 2)
      .attr("y", height + margin.top + 35)
      .attr("text-anchor", "middle")
      .attr("font-size", "12px")
      .text(segmentB);
  
    svg.append("text")
      .attr("transform", `rotate(-90)`)
      .attr("x", -margin.top - height / 2)
      .attr("y", 12)
      .attr("text-anchor", "middle")
      .attr("font-size", "12px")
      .text(segmentA);
  
    // Zoom behavior
    zoom = d3.zoom()
      .scaleExtent([1, 10])
      .translateExtent([[-100, -100], [width + 100, height + 100]])
      .on("zoom", zoomed);
  
    svg.call(zoom);
  
    function zoomed({ transform }) {
      // Zoom heatmap
      heatmap.attr("transform", transform);
  
      // Zoom axes
      const zx = transform.rescaleX(xScale);
      const zy = transform.rescaleY(yScale);
  
      gx.call(xAxis.scale(zx));
      gy.call(yAxis.scale(zy));
    }

  };

  const renderHeatmap = () => {
    // make request to backend
    const host = process.env.NEXT_PUBLIC_DJANGO_HOST;
    fetch(`${host}/api/timge/heatmap/`, {
      method: "POST",
      headers: {
        "Content-Type": "application/json"
      },
      body: JSON.stringify({
        uuid: "1234",
        file_name: track,
        genome_path: reference,
        resolution: resolution,
        segment_1: segmentA,
        segment_2: segmentB
      })
    })
      .then(response => response.json())
      .then(data => {
        if (data.status === "success" && data.matrix) {
          drawHeatmap(data.matrix);
        } else {
          console.error("Failed to generate heatmap", data.message);
        }
      }
      )
  }

  return (
    <ParentView viewConfig={props.viewConfig}>
        {
            ! (reference && track) ? (
        <Box
            sx={{
                display: "flex",
                flexDirection: "row",
                justifyContent: "center",
                alignItems: "center",
                width: "100%",
            }}
        >
            <Typography fontSize="md" margin={"0 10px"}>
            Select reference
            </Typography>
            <Select
            defaultValue={reference}
            placeholder="Select reference genome"
            onChange={(e, value) => setReference(value)}
            sx={{
                margin: "0 10px",}}
            >
            {
                props.trackFiles.map((file, index) => (
                    <Option key={index} value={file.name}>
                    {file.name}
                    </Option>
                ))}
            </Select>
            <Typography fontSize="md" margin={"0 10px"}>
            Select track
            </Typography>
            <Select
            defaultValue={track}
            placeholder="Select track"
            onChange={(e, value) => setTrack(value)}
            sx={{
                margin: "0 10px",}}
            >
            {
                props.trackFiles.map((file, index) => (
                    <Option key={index} value={file.name}>
                    {file.name}
                    </Option>
                ))}
            </Select>
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
                    onChange={(e, value) => setSegmentA(value)}
                    defaultValue={segmentA}
                    placeholder="Select segment A"
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
                    onChange={(e, value) => setSegmentB(value)}
                    defaultValue={segmentB}
                    placeholder="Select segment B"
                    >
                    {availableSegments.map((segment, index) => (
                        <Option key={index} value={segment}>
                        {segment}
                        </Option>
                    ))}
                    </Select>

                    <Typography fontSize="md">
                    Resolution:
                    </Typography>
                    <Select
                    defaultValue={resolution}
                    placeholder="Select resolution"
                    onChange={(e, value) => setResolution(value)}
                    >
                    <Option value={5}>5</Option>
                    <Option value={25}>25</Option>
                    <Option value={50}>50</Option>
                    <Option value={100}>100</Option>
                    <Option value={200}>200</Option>
                    <Option value={500}>500</Option>
                    <Option value={1000}>1000</Option>
                    </Select>

                    <Button variant="solid" color="primary" onClick={renderHeatmap}>
                    Render
                    </Button>
                    <Button
                    variant="outlined"
                    color="neutral"
                    onClick={() => {
                        const svg = d3.select(heatmapRef.current);
                        svg.transition()
                        .duration(500)
                        .call(zoom.transform, d3.zoomIdentity); // Reset zoom!
                    }}
                    >
                    Reset Zoom
                </Button>
                </Box>
            </Card>
            <Box
            sx={{
                margin: "1em",
            }}
            >
                <svg ref={heatmapRef}></svg>
            </Box>
        </Box>
            )
        }
    </ParentView>

  );
};

export default MapView;
