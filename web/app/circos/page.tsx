"use client"

import { useEffect, useRef, useState } from 'react';
import * as d3 from 'd3';
import { Assembly } from '../types/genomes';
import FileUpload from '../../components/FileUpload';
import Form from '../../components/Form';

const CanvasPage = () => {
  const canvasRef = useRef<HTMLDivElement>(null);
  const [segments, setSegments] = useState<Assembly[]>([]);
  const [segmentAnglePadding, setSegmentAnglePadding] = useState(0.02);
  const [axisLabelFontSize, setAxisLabelFontSize] = useState(6);
  const [showAxis, setShowAxis] = useState(true);
  const [segmentInnerRadius, setSegmentInnerRadius] = useState(200);
  const [segmentOuterRadius, setSegmentOuterRadius] = useState(250);
  const [showAxisLabels, setShowAxisLabels] = useState(true);
  const [canvasWidth, setCanvasWidth] = useState(800);
  const [canvasHeight, setCanvasHeight] = useState(600);
  const [segmentGridPadding, setSegmentGridPadding] = useState(3);
  const [tickLength, setTickLength] = useState(3);
  const [tickTextPadding, setTickTextPadding] = useState(3);

  const color_palette = d3.scaleSequential(d3.interpolateViridis);
  const tick_precision = 200000;


  const handleFileUpload = (file: File) => {
    const reader = new FileReader();
    reader.onload = (event) => {
      try {
        const json = JSON.parse(event.target?.result as string);
        setSegments(json);
      } catch (error) {
        console.error('Invalid JSON file:', error);
        alert('Uploaded file is not valid JSON.');
      }
    };
    reader.readAsText(file);
  };

  const updateAssemblyConfig = (
    segmentPadding: number,
    axisLabelFontSize: number,
    showAxis: boolean, 
    segmentInnerRadius: number, 
    segmentOuterRadius: number,
    segmentGridPadding: number,
    showAxisLabels: boolean,
    canvasWidth: number,
    canvasHeight: number,
    tickLength: number,
    tickTextPadding: number
  ) => {
    setSegmentAnglePadding(segmentPadding);
    setAxisLabelFontSize(axisLabelFontSize);
    setShowAxis(showAxis);
    setSegmentInnerRadius(segmentInnerRadius);
    setSegmentOuterRadius(segmentOuterRadius);
    setSegmentGridPadding(segmentGridPadding);
    setShowAxisLabels(showAxisLabels);
    setCanvasWidth(canvasWidth);
    setCanvasHeight(canvasHeight);
    setTickLength(tickLength);
    setTickTextPadding(tickTextPadding);
  };

  useEffect(() => {
    if (canvasRef.current && segments.length > 0) {
      const total_length = segments.reduce((acc, segment) => acc + (segment.end - segment.start), 0);
      const total_angle_padding = 2 * segmentAnglePadding * segments.length;
      const total_available_angle = 2 * Math.PI - total_angle_padding;

      d3.select(canvasRef.current).select("svg").remove();

      const svg = d3
        .select(canvasRef.current)
        .append("svg")
        .attr("viewBox", [0, 0, canvasWidth, canvasHeight])
        .attr("width", canvasWidth)
        .attr("height", canvasHeight)
        .style("border", "1px solid black");

      let last_angle = 0;

      segments.forEach((segment, i) => {
        const startAngle = last_angle + segmentAnglePadding;
        const endAngle = last_angle + segmentAnglePadding + ((segment.end - segment.start) / total_length) * total_available_angle;

        const arc = d3.arc()
          .innerRadius(segmentInnerRadius)
          .outerRadius(segmentOuterRadius)
          .startAngle(startAngle)
          .endAngle(endAngle);



        svg.append("path")
          .attr("d", arc)
          .attr("transform", `translate(${canvasWidth / 2}, ${canvasHeight / 2})`)
          .attr("fill", color_palette(0.55 * i/segments.length))
          .on("mouseover", function (this: SVGPathElement) {
            d3.select(this).attr("filter", "brightness(0.9)");
          })
          .on("mouseout", function (this: SVGPathElement) {
            d3.select(this).attr("filter", "brightness(1)");
          });
        
          svg.append("text")
          .attr("x", canvasWidth / 2)
          .attr("y", canvasHeight / 2)
          .attr("text-anchor", "middle")
          .attr("dominant-baseline", "central")
          .attr("font-size", 12)
          .attr("transform", `translate(${(segmentInnerRadius + segmentOuterRadius) / 2 * Math.sin((startAngle + endAngle) / 2)}, ${-(segmentInnerRadius + segmentOuterRadius) / 2 * Math.cos((startAngle + endAngle) / 2)})`)
          .attr("fill", "white")
          .text(segment.chromosome);

        if (showAxis) {
            const arc_grid = d3.arc()
                .innerRadius(segmentOuterRadius + segmentGridPadding)
                .outerRadius(segmentOuterRadius + segmentGridPadding)
                .startAngle(startAngle)
                .endAngle(endAngle);
        
            svg.append("path")
                .attr("d", arc_grid)
                .attr("transform", `translate(${canvasWidth / 2}, ${canvasHeight / 2})`)
                .attr("stroke", "black")
                .attr("stroke-width", 1);
            
            const numTicks = Math.floor((segment.end - segment.start) / tick_precision);
            const scale = d3.scaleLinear()
            .domain([segment.start, segment.end])
            .range([startAngle, endAngle]);

            const ticks = scale.ticks(numTicks);

            ticks.forEach((tick: any) => {
                const angle = scale(tick);
                svg
                    .append("line")
                    .attr("x1", (canvasWidth / 2) + (segmentOuterRadius + segmentGridPadding) * Math.sin(angle))
                    .attr("y1", (canvasHeight / 2) - (segmentOuterRadius + segmentGridPadding) * Math.cos(angle))
                    .attr("x2", (canvasWidth / 2) + (segmentOuterRadius + segmentGridPadding + tickLength) * Math.sin(angle))
                    .attr("y2", (canvasHeight / 2) - (segmentOuterRadius + segmentGridPadding + tickLength) * Math.cos(angle))
                    .attr("stroke", "black")
                    .attr("stroke-width", 1);
                
                const x = (canvasWidth / 2) + (segmentOuterRadius + segmentGridPadding + tickLength + tickTextPadding) * Math.sin(angle);
                const y = (canvasHeight / 2) - (segmentOuterRadius + segmentGridPadding + tickLength + tickTextPadding) * Math.cos(angle);
                const rotation = ((angle * 180) / Math.PI) - 90;
                const tick_text = tick == 0 ? "" : tick.toLocaleString();
                svg
                    .append("text")
                    .attr("x", x)
                    .attr("y", y)
                    .attr("text-anchor", "start")
                    .attr("dominant-baseline", "central")
                    .attr("font-size", axisLabelFontSize)
                    .attr("transform", `rotate(${rotation}, ${x}, ${y})`)
                    .text(tick_text);
            });
        }

        last_angle = endAngle + segmentAnglePadding;
      });
    }
  }, [
    segments, 
    segmentAnglePadding, 
    axisLabelFontSize, 
    showAxis, 
    segmentInnerRadius, 
    segmentOuterRadius, 
    segmentGridPadding, 
    showAxisLabels, 
    canvasWidth, 
    canvasHeight,
    tickLength,
    tickTextPadding
]);

  return (
    <div>
      <h1>D3 Canvas Page</h1>
      <FileUpload onFileUpload={handleFileUpload} />
      <div ref={canvasRef}></div>
      {segments.length > 0 ? <Form onUpdate={updateAssemblyConfig}/> : null}
    </div>
  );
};

export default CanvasPage;
