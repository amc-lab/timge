import React, { useRef, useEffect, useState } from "react";
import * as d3 from "d3";
import { View } from "@/store/features/views/types";

interface CanvasHeatmapProps {
  matrix: number[][];
  segmentA: string;
  segmentB: string;
  title?: string;
  resolution: number;
  toggleColourScheme: boolean;
  showGridlines: boolean;
  isMinimised: boolean;
  setCanvasRef?: (el: HTMLCanvasElement | null) => void;
  setZoomRef?: (zoom: d3.ZoomBehavior<Element, unknown>, svgEl: SVGSVGElement | null) => void;
  viewConfig?: View;
  onLocusChange?: (xRange: [number, number], yRange: [number, number]) => void;
  zoomToLocusRef?: React.MutableRefObject<((x0: number, x1: number, y0: number, y1: number) => void) | null>;
}

const CanvasHeatmap = ({
  matrix,
  title,
  segmentA,
  segmentB,
  resolution,
  toggleColourScheme,
  showGridlines,
  isMinimised,
  setCanvasRef,
  setZoomRef,
  viewConfig,
  onLocusChange,
  zoomToLocusRef,
}: CanvasHeatmapProps) => {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const svgRef = useRef<SVGSVGElement>(null);
  const gxRef = useRef<SVGGElement>(null);
  const gyRef = useRef<SVGGElement>(null);
  const xScaleRef = useRef<d3.ScaleLinear<number, number> | null>(null);
  const yScaleRef = useRef<d3.ScaleLinear<number, number> | null>(null);
  const transformRef = useRef(d3.zoomIdentity);

  const margin = { top: 40, right: 70, bottom: 70, left: 70 };
  const fixedWidth = isMinimised ? 400 : 800;

  useEffect(() => {
    if (!matrix || matrix.length === 0) return;

    if (setCanvasRef) {
      setCanvasRef(canvasRef.current);
    }

    const canvas = canvasRef.current!;
    const ctx = canvas.getContext("2d")!;
    ctx.imageSmoothingEnabled = false;
    const svg = d3.select(svgRef.current)
      .style("font-family", "Arial");
    const gx = d3.select(gxRef.current);
    const gy = d3.select(gyRef.current);
    svg.selectAll("defs").remove();

    const numRows = matrix.length;
    const numCols = matrix[0].length;

    // const width = canvasRef.current.clientWidth;
    const width = fixedWidth;<s></s>
    // const cellSize = fixedWidth / numRows;
    const cellSize = width / numRows; 
    const height = cellSize * numCols;  

    canvas.width = width;
    canvas.height = height;
    canvas.style.width = `${width}px`;
    canvas.style.height = `${height}px`;

    svg.attr("width", width + margin.left + margin.right)
       .attr("height", height + margin.top + margin.bottom);

    const xScale = d3.scaleLinear().domain([0, numRows]).range([0, width]);
    const yScale = d3.scaleLinear().domain([0, numCols]).range([height, 0]);
    xScaleRef.current = xScale;
    yScaleRef.current = yScale;

      // remove old title
    svg.selectAll(".chart-title").remove();

    // draw new title
    if (title) {
      svg.append("text")
        .attr("class", "chart-title")
        // center across the full SVG width
        .attr("x", (width + margin.left + margin.right) / 2)
        // place halfway through the top margin
        .attr("y", margin.top / 2)
        .attr("text-anchor", "middle")
        .attr("font-size", "16px")
        .attr("font-weight", "bold")
        .text(title);
    }

    gx.attr("transform", `translate(${margin.left},${margin.top + height})`)
      .call(d3.axisBottom(xScale).ticks(10).tickFormat(d => `${+d * resolution}`));

    gy.attr("transform", `translate(${margin.left},${margin.top})`)
      .call(d3.axisLeft(yScale).ticks(10).tickFormat(d => `${+d * resolution}`));

    svg.selectAll(".axis-label").remove();

    svg.append("text")
      .attr("class", "axis-label")
      .attr("x", margin.left + width / 2)
      .attr("y", height + margin.top + 50)
      .attr("text-anchor", "middle")
      .attr("font-size", "12px")
      .attr("font-family", "Arial")
      .text(segmentA);

    svg.append("text")
      .attr("class", "axis-label")
      .attr("transform", `rotate(-90)`)
      .attr("x", -margin.top - height / 2)
      .attr("y", 15)
      .attr("text-anchor", "middle")
      .attr("font-size", "12px")
      .attr("font-family", "Arial")
      .text(segmentB);

    const colorScale = d3.scaleSequential()
      .domain([0, d3.max(matrix.flat())!])
      .interpolator(toggleColourScheme ? d3.interpolateOrRd : d3.interpolateGnBu);

    const legendHeight = 200;
    const legendWidth = 15;

    const legendScale = d3.scaleLinear()
      .domain(colorScale.domain())
      .range([legendHeight, 0]);

      svg.selectAll(".heatmap-legend").remove();
      svg.select("defs").remove();
      
      const defs = svg.append("defs");
      const gradientId = "legend-gradient";
      
      const linearGradient = defs.append("linearGradient")
        .attr("id", gradientId)
        .attr("x1", "0%").attr("y1", "100%")
        .attr("x2", "0%").attr("y2", "0%");
      
      const numStops = 10;
      const step = 1 / (numStops - 1);
      d3.range(numStops).forEach(i => {
        linearGradient.append("stop")
          .attr("offset", `${i * step * 100}%`)
          .attr("stop-color", colorScale(colorScale.domain()[0] + i * step * (colorScale.domain()[1] - colorScale.domain()[0])));
      });
      
      const legendGroup = svg.append("g")
        .attr("class", "heatmap-legend")
        .attr("transform", `translate(${margin.left + width + 30}, ${margin.top})`);
      
      legendGroup.append("rect")
        .attr("width", 15)
        .attr("height", 200)
        .style("fill", `url(#${gradientId})`);
      
      legendGroup.append("g")
        .attr("transform", `translate(15, 0)`)
        .call(d3.axisRight(legendScale).ticks(5));

    const zoom = d3.zoom()
      .scaleExtent([1, 10])
      .translateExtent([[0, 0], [width, height]])
      .extent([[0, 0], [width, height]])
      .on("zoom", ({ transform }) => {
        transformRef.current = transform;

        const zx = transform.rescaleX(xScale);
        const zy = transform.rescaleY(yScale);

        if (onLocusChange) {
          const zx = transform.rescaleX(xScale);
          const zy = transform.rescaleY(yScale);
          const xDomain = zx.domain().map(d => Math.round(d * resolution)) as [number, number];
          const yDomain = zy.domain().map(d => Math.round(d * resolution)) as [number, number];
          onLocusChange(xDomain, yDomain);
        }


        gx.call(d3.axisBottom(zx).ticks(10).tickFormat(d => `${Math.round(+d * resolution)}`));
        gy.call(d3.axisLeft(zy).ticks(10).tickFormat(d => `${Math.round(+d * resolution)}`));

        draw();
      });
    
    if (setZoomRef) {
      setZoomRef(zoom, svgRef.current);
    }

    if (zoomToLocusRef) {
      zoomToLocusRef.current = (x0, x1, y0, y1) => {
        const scaleX = fixedWidth / ((x1 - x0) / resolution);
        const scaleY = (fixedWidth * matrix[0].length / matrix.length) / ((y1 - y0) / resolution);
        const scale = Math.min(scaleX, scaleY);
        const tx = -x0 / resolution * scale;
        const ty = -y0 / resolution * scale;
        const transform = d3.zoomIdentity.translate(tx, ty).scale(scale);
        d3.select(svgRef.current).call(zoom.current!.transform, transform);
      };
    }

    svg.call(zoom as any);

    const draw = () => {
      ctx.clearRect(0, 0, canvas.width, canvas.height);
      ctx.save();
      const t = transformRef.current;
      ctx.translate(t.x + width * t.k, t.y);
      ctx.scale(-t.k, t.k);
    
      for (let i = 0; i < numRows; i++) {
        for (let j = 0; j < numCols; j++) {
          ctx.fillStyle = colorScale(matrix[i][j])!;
          const x = Math.floor(i * cellSize);
          const y = Math.floor(j * cellSize);
          const size = Math.ceil(cellSize);
    
          ctx.fillRect(x, y, size, size);
    
          // Only stroke if explicitly enabled
          if (showGridlines) {
            ctx.strokeStyle = "#000"; // or any color
            ctx.strokeRect(x, y, size, size);
          }
        }
      }
    
      ctx.restore();
    };
    

    draw();
  }, [canvasRef, matrix, toggleColourScheme, showGridlines, resolution, isMinimised, segmentA, segmentB]);

  return (
    <div
      style={{
        position: "relative",
        width: `${fixedWidth + margin.left + margin.right}px`,
        height: "auto",
      }}
    >
      <canvas
        ref={canvasRef}
        style={{
          position: "absolute",
          top: `${margin.top}px`,
          left: `${margin.left}px`,
          zIndex: 0,
        }}
      />
      <svg
        ref={svgRef}
        style={{
          position: "relative",
          zIndex: 1,
        }}
      >
        <g ref={gxRef} className="x-axis" />
        <g ref={gyRef} className="y-axis" />
      </svg>
    </div>
  );
};

export default CanvasHeatmap;
