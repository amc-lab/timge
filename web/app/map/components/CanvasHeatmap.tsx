import React, { useRef, useEffect, useState } from "react";
import * as d3 from "d3";

interface CanvasHeatmapProps {
  matrix: number[][];
  segmentA: string;
  segmentB: string;
  resolution: number;
  toggleColourScheme: boolean;
  showGridlines: boolean;
  isMinimised: boolean;
}

const CanvasHeatmap = ({
  matrix,
  segmentA,
  segmentB,
  resolution,
  toggleColourScheme,
  showGridlines,
  isMinimised,
}: CanvasHeatmapProps) => {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const svgRef = useRef<SVGSVGElement>(null);
  const gxRef = useRef<SVGGElement>(null);
  const gyRef = useRef<SVGGElement>(null);
  const xScaleRef = useRef<d3.ScaleLinear<number, number> | null>(null);
  const yScaleRef = useRef<d3.ScaleLinear<number, number> | null>(null);
  const transformRef = useRef(d3.zoomIdentity);

  const margin = { top: 40, right: 70, bottom: 70, left: 70 };
  const fixedWidth = isMinimised ? 500 : 800;

  useEffect(() => {
    if (!matrix || matrix.length === 0) return;

    const canvas = canvasRef.current!;
    const ctx = canvas.getContext("2d")!;
    const svg = d3.select(svgRef.current);
    const gx = d3.select(gxRef.current);
    const gy = d3.select(gyRef.current);
    svg.selectAll("defs").remove(); // clear any old defs

    const numRows = matrix.length;
    const numCols = matrix[0].length;

    const cellSize = fixedWidth / numRows;
    const width = fixedWidth;
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
      .text(segmentA);

    svg.append("text")
      .attr("class", "axis-label")
      .attr("transform", `rotate(-90)`)
      .attr("x", -margin.top - height / 2)
      .attr("y", 15)
      .attr("text-anchor", "middle")
      .attr("font-size", "12px")
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

        gx.call(d3.axisBottom(zx).ticks(10).tickFormat(d => `${Math.round(+d * resolution)}`));
        gy.call(d3.axisLeft(zy).ticks(10).tickFormat(d => `${Math.round(+d * resolution)}`));

        draw();
      });

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
          ctx.fillRect(i * cellSize, j * cellSize, cellSize, cellSize);
          if (showGridlines) {
            ctx.strokeStyle = "#ccc";
            ctx.strokeRect(i * cellSize, j * cellSize, cellSize, cellSize);
          }
        }
      }

      ctx.restore();
    };

    draw();
  }, [matrix, toggleColourScheme, showGridlines, resolution, isMinimised, segmentA, segmentB]);

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
