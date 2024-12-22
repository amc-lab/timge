"use client"
import { useEffect, useRef } from 'react';
import * as d3 from 'd3';

const CanvasPage = () => {
  const canvasRef = useRef<HTMLDivElement>(null);

  const segments = [
    {
      "chromosome": "os1",
      "id": 1,
      "start": 0,
      "end": 640851,
      "color": "yellow"
    },
    {
      "chromosome": "os2",
      "id": 2,
      "start": 0,
      "end": 947102,
      "color": "gpos25"
    },
    {
      "chromosome": "os3",
      "id": 3,
      "start": 0,
      "end": 1067971,
      "color": "black"
    },
    {
      "chromosome": "os4",
      "id": 4,
      "start": 0,
      "end": 1200490,
      "color": "blue"
    },
    {
      "chromosome": "os5",
      "id": 5,
      "start": 0,
      "end": 1343557,
      "color": "gpos75"
    },
    {
      "chromosome": "os6",
      "id": 6,
      "start": 0,
      "end": 1418242,
      "color": "spectral-5-div-4"
    },
    {
      "chromosome": "os7",
      "id": 7,
      "start": 0,
      "end": 1445207,
      "color": "spectral-5-div-4"
    },
    {
      "chromosome": "os8",
      "id": 8,
      "start": 0,
      "end": 1472805,
      "color": "gpos100"
    },
    {
      "chromosome": "os9",
      "id": 9,
      "start": 0,
      "end": 1541735,
      "color": "stalk"
    },
    {
      "chromosome": "os10",
      "id": 10,
      "start": 0,
      "end": 1687656,
      "color": "acen"
    },
    {
      "chromosome": "os11",
      "id": 11,
      "start": 0,
      "end": 2038340,
      "color": "green"
    },
    {
      "chromosome": "os12",
      "id": 12,
      "start": 0,
      "end": 2271494,
      "color": "orange"
    },
    {
      "chromosome": "os13",
      "id": 13,
      "start": 0,
      "end": 2925236,
      "color": "blue"
    },
    {
      "chromosome": "os14",
      "id": 14,
      "start": 0,
      "end": 3291936,
      "color": "spectral-5-div-1"
    }
  ];

  const total_length = segments.reduce((acc, segment) => acc + segment.end, 0);
  const color_palette = d3.scaleOrdinal(d3.schemeTableau10);
  const angle_padding = 0.01;
  const total_angle_padding = 2 * angle_padding * segments.length;
  const total_available_angle = 2 * Math.PI - total_angle_padding;

  const tick_precision = 10000;
  const total_ticks = Math.floor(total_length / tick_precision);
  const tick_angle = total_available_angle / total_ticks;

  console.log(total_length, total_available_angle);

  useEffect(() => {
    if (canvasRef.current) {
      d3.select(canvasRef.current).select("svg").remove();

      const svg = d3
        .select(canvasRef.current)
        .append("svg")
        .attr("width", 800)
        .attr("height", 600)
        .style("border", "1px solid black");

      let last_angle = 0;

      segments.forEach((segment, i) => {
        const startAngle = last_angle + angle_padding;
        const endAngle = last_angle + angle_padding + ((segment.end - segment.start) / total_length) * total_available_angle;

        const arc = d3.arc()
          .innerRadius(200)
          .outerRadius(250)
          .startAngle(startAngle)
          .endAngle(endAngle);

        const arc_grid = d3.arc()
            .innerRadius(253)
            .outerRadius(254)
            .startAngle(startAngle)
            .endAngle(endAngle);

        last_angle = endAngle + angle_padding;

        svg.append("path")
          .attr("d", arc)
          .attr("transform", "translate(400, 300)")
          .attr("fill", color_palette(i))
          .on("mouseover", function (this: SVGPathElement) {
            d3.select(this).attr("filter", "brightness(0.9)");
          })
          .on("mouseout", function (this: SVGPathElement) {
            d3.select(this).attr("filter", "brightness(1)");
          });

        svg.append("path")
            .attr("d", arc_grid)
            .attr("transform", "translate(400, 300)")
            .attr("fill", "black");

      });
    }
  }, []);

  return (
    <div>
      <h1>D3 Canvas Page</h1>
      <div ref={canvasRef}></div>
    </div>
  );
};

export default CanvasPage;
