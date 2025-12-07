import React, { useRef, useEffect } from "react";
import * as d3 from "d3";
import { buildTileCacheKey, TILE_SIZE_BINS, tileCache } from "../../map/TileCache";
import { Segment } from '../MapView';

interface CanvasHeatmapProps {
  segmentA: Segment | null;
  segmentB: Segment | null;
  title?: string;
  trackName: string;
  resolution: number;
  tileCacheVersion: number;
  toggleColourScheme: boolean;
  isMinimised: boolean;
  setCanvasRef?: (el: HTMLCanvasElement | null) => void;
  setZoomRef?: (zoom: d3.ZoomBehavior<Element, unknown>, svgEl: SVGSVGElement | null) => void;
  onLocusChange?: (xRange: [number, number], yRange: [number, number]) => void;
  onZoomEnd?: (xRange: [number, number], yRange: [number, number]) => void;
  zoomToLocusRef?: React.MutableRefObject<
    | ((x0: number, x1: number, y0: number, y1: number) => void)
    | null
  >;
}

const TILE_BINS = TILE_SIZE_BINS;

const CanvasHeatmap = ({
  segmentA,
  segmentB,
  title,
  trackName,
  resolution,
  tileCacheVersion,
  toggleColourScheme,
  isMinimised,
  setCanvasRef,
  setZoomRef,
  onLocusChange,
  onZoomEnd,
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

  // ------------------------------------------------------------
  //  MAIN EFFECT – initialize canvas, axes, zoom and drawing
  // ------------------------------------------------------------
  useEffect(() => {
    if (!segmentA?.id || !segmentB?.id || !resolution) {
      const canvas = canvasRef.current;
      const ctx = canvas?.getContext("2d");
      if (ctx && canvas) {
        ctx.clearRect(0, 0, canvas.width, canvas.height);
      }
      return;
    }

    console.log("CanvasHeatmap render effect", {
      segmentA,
      segmentB,
      resolution,
      toggleColourScheme,
      isMinimised,
    });
    if (setCanvasRef) setCanvasRef(canvasRef.current);

    const canvas = canvasRef.current!;
    const ctx = canvas.getContext("2d")!;
    ctx.imageSmoothingEnabled = false;

    const svg = d3.select(svgRef.current);
    const gx = d3.select(gxRef.current);
    const gy = d3.select(gyRef.current);

    const width = fixedWidth - 40;
    // const height = width; // square canvas (same as original behaviour)
    const height = width * (segmentB.length / segmentA.length)

    canvas.width = width;
    canvas.height = height;
    canvas.style.width = `${width}px`;
    canvas.style.height = `${height}px`;

    svg.attr("width", width + margin.left + margin.right + 40)
       .attr("height", height + margin.top + margin.bottom);

    // ----------------------------------------------
    // AXIS DOMAIN IS ALWAYS IN BINS:
    //   x = [0 .. total_bins_x]
    //   y = [0 .. total_bins_y]
    // BACKEND defines tile coverage, not here.
    // We treat the entire canvas as [0..N] where N is
    // whatever tiles cover (max tile index * TILE_BINS)
    // ----------------------------------------------
    const keyPrefix = `${trackName}|${segmentA.id}|${segmentB.id}|${resolution}|`;
    const relevantKeys = tileCache.keys().filter((k) => k.startsWith(keyPrefix));
    console.log("Tile cache keys for segments:", relevantKeys);

    const totalBinsX = Math.ceil(segmentA.length / resolution);
    const totalBinsY = Math.ceil(segmentB.length / resolution);

    const xScale = d3.scaleLinear().domain([0, totalBinsX]).range([0, width]);
    const yScale = d3.scaleLinear().domain([0, totalBinsY]).range([height, 0]);

    xScaleRef.current = xScale;
    yScaleRef.current = yScale;

    // ------------------------------------------
    // TITLE
    // ------------------------------------------
    svg.selectAll(".chart-title").remove();
    if (title) {
      svg.append("text")
        .attr("class", "chart-title")
        .attr("x", (width + margin.left + margin.right) / 2)
        .attr("y", margin.top / 2)
        .attr("text-anchor", "middle")
        .attr("font-size", "16px")
        .attr("font-weight", "bold")
        .text(title);
    }

    // ------------------------------------------
    // AXES
    // ------------------------------------------
    gx.attr("transform", `translate(${margin.left},${margin.top + height})`)
      .call(d3.axisBottom(xScale).ticks(10).tickFormat(d => `${d * resolution}`));

    gy.attr("transform", `translate(${margin.left},${margin.top})`)
      .call(d3.axisLeft(yScale).ticks(10).tickFormat(d => `${d * resolution}`));

    // Axis labels
    svg.selectAll(".axis-label").remove();
    svg.append("text")
      .attr("class", "axis-label")
      .attr("x", margin.left + width / 2)
      .attr("y", height + margin.top + 50)
      .attr("text-anchor", "middle")
      .text(segmentA.id);

    svg.append("text")
      .attr("class", "axis-label")
      .attr("transform", "rotate(-90)")
      .attr("x", -margin.top - height / 2)
      .attr("y", 15)
      .attr("text-anchor", "middle")
      .text(segmentB.id);

    if (onLocusChange) {
      const xBP = xScale.domain().map((d) => Math.round(d * resolution)) as [number, number];
      const yBP = yScale.domain().map((d) => Math.round(d * resolution)) as [number, number];
      onLocusChange(xBP, yBP);
    }

    // ------------------------------------------
    // COLOUR SCALE (based on all tile values)
    // ------------------------------------------
    const allValues: number[] = [];
    for (const key of relevantKeys) {
      const tile = tileCache.get(key);
      if (tile) allValues.push(...(tile as number[][]).flat());
    }

    const vmax = d3.max(allValues) ?? 1;
    const colorScale = d3.scaleSequential()
      .domain([0, vmax])
      .interpolator(toggleColourScheme ? d3.interpolateOrRd : d3.interpolateGnBu);

    // ------------------------------------------
    // LEGEND
    // ------------------------------------------
    svg.selectAll(".heatmap-legend").remove();
    svg.select("defs").remove();

    const defs = svg.append("defs");
    const gradientId = "legend-gradient";

    const linearGradient = defs.append("linearGradient")
      .attr("id", gradientId)
      .attr("x1", "0%")
      .attr("y1", "100%")
      .attr("x2", "0%")
      .attr("y2", "0%");

    const numStops = 10;
    const step = 1 / (numStops - 1);
    for (let i = 0; i < numStops; i++) {
      linearGradient.append("stop")
        .attr("offset", `${i * step * 100}%`)
        .attr("stop-color",
          colorScale(i * step * vmax)
        );
    }

    const legendGroup = svg.append("g")
      .attr("class", "heatmap-legend")
      .attr("transform", `translate(${margin.left + width + 30}, ${margin.top})`);

    const legendHeight = 200;
    legendGroup.append("rect")
      .attr("width", 15)
      .attr("height", legendHeight)
      .style("fill", `url(#${gradientId})`);

    const legendScale = d3.scaleLinear()
      .domain(colorScale.domain())
      .range([legendHeight, 0]);

    legendGroup.append("g")
      .attr("transform", "translate(15,0)")
      .call(d3.axisRight(legendScale).ticks(5));

    // ------------------------------------------
    // ZOOM
    // ------------------------------------------
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

        if (onLocusChange) {
          const xBP = zx.domain().map(d => Math.round(d * resolution)) as [number, number];
          const yBP = zy.domain().map(d => Math.round(d * resolution)) as [number, number];
          onLocusChange(xBP, yBP);
        }

        draw();
      })
      .on("end", ({ transform }) => {
        if (!onZoomEnd) return;
        const zx = transform.rescaleX(xScale);
        const zy = transform.rescaleY(yScale);
        const xBP = zx.domain().map((d) => Math.round(d * resolution)) as [number, number];
        const yBP = zy.domain().map((d) => Math.round(d * resolution)) as [number, number];
        onZoomEnd(xBP, yBP);
      });

    svg.call(zoom as any);

    if (setZoomRef) setZoomRef(zoom, svgRef.current);

    // allow external zoomToLocus()
    if (zoomToLocusRef) {
      zoomToLocusRef.current = (x0, x1, y0, y1) => {
        const scaleX = width / ((x1 - x0) / resolution);
        const scaleY = height / ((y1 - y0) / resolution);
        const scale = Math.min(scaleX, scaleY);

        const tx = -x0 / resolution * scale;
        const ty = -y0 / resolution * scale;

        const transform = d3.zoomIdentity.translate(tx, ty).scale(scale);
        d3.select(svgRef.current).call(zoom.transform, transform);
      };
    }

    // ------------------------------------------
    // DRAW FUNCTION — TILE RENDERING
    // ------------------------------------------
    function draw() {
      ctx.clearRect(0, 0, width, height);
      ctx.save();

      const T = transformRef.current;

      ctx.translate(T.x, T.y);
      ctx.scale(T.k, T.k);

      const invK = 1 / T.k;
      const viewX0 = xScale.invert(-T.x * invK);
      const viewX1 = xScale.invert((width - T.x) * invK);

      const viewY0 = yScale.invert((height - T.y) * invK);
      const viewY1 = yScale.invert(-T.y * invK);

      const tileX0 = Math.floor(viewX0 / TILE_BINS);
      const tileX1 = Math.floor(viewX1 / TILE_BINS);
      const tileY0 = Math.floor(viewY0 / TILE_BINS);
      const tileY1 = Math.floor(viewY1 / TILE_BINS);

      for (let tx = tileX0; tx <= tileX1; tx++) {
        for (let ty = tileY0; ty <= tileY1; ty++) {
          const key = buildTileCacheKey({
            trackName,
            segmentAId: segmentA.id,
            segmentBId: segmentB.id,
            resolution,
            tileX: tx,
            tileY: ty,
          });
          const tile = tileCache.get(key);
          if (!tile) continue;

          for (let i = 0; i < TILE_BINS; i++) {
            for (let j = 0; j < TILE_BINS; j++) {
              const binValue = tile[i][j];
              ctx.fillStyle = colorScale(binValue)!;

              const x = xScale(tx * TILE_BINS + i);
              const y = yScale(ty * TILE_BINS + j + 1);

              const x2 = xScale(tx * TILE_BINS + i + 1);
              const y2 = yScale(ty * TILE_BINS + j);

              const w = x2 - x;
              const h = y - y2;

              ctx.fillRect(x, y2, w, h);
            }
          }
        }
      }

      ctx.restore();
    }

    draw();
  }, [
    tileCacheVersion,
    trackName,
    resolution,
    toggleColourScheme,
    isMinimised,
    segmentA,
    segmentB,
  ]);

  // ------------------------------------------------------------
  //  RENDER WRAPPER CONTAINERS
  // ------------------------------------------------------------
  return (
    <div
      style={{
        position: "relative",
        width: `${fixedWidth + margin.left + margin.right}px`,
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
