import { useEffect } from "react";
import * as d3 from "d3";
import { BarData, BarConfig, GlobalConfig } from "../../types/genomes";

interface BarProps {
  data: {
    bars: Array<BarData>;
    globalConfig: GlobalConfig;
    divRef: any;
  };
  config: BarConfig;
  segments: Array<any>;
  idx: number;
}

const Bar = ({ data, config, segments, idx }: BarProps) => {
  const canvasRef = data.divRef;

  useEffect(() => {
    let svg = d3.select(canvasRef.current).select("svg");

    if (svg.empty()) {
      svg = d3
        .select(canvasRef.current)
        .append("svg")
        .attr("width", data.globalConfig.canvasWidth)
        .attr("height", data.globalConfig.canvasHeight);
    }

    const uniqueGroupClass = `group-${idx}`;

    svg.selectAll(`g.${uniqueGroupClass}`).remove();

    const group = svg
      .append("g")
      .attr("class", uniqueGroupClass)
      .attr(
        "transform",
        `translate(${data.globalConfig.canvasWidth / 2}, ${data.globalConfig.canvasHeight / 2})`,
      );

    const innerRadius = config.innerRadius;
    const outerRadius = config.innerRadius + 0.8 * config.trackWidth;

    // Scale for bar height
    const barScale = d3
      .scaleLinear()
      .domain([0, d3.max(data.bars, (d) => d.value) || 100])
      .range([innerRadius, outerRadius]);

    // Draw bars within segment spaces
    group
      .selectAll("path.bar")
      .data(data.bars)
      .enter()
      .append("path")
      .attr("class", "bar")
      .attr("d", (d) => {
        const segment = segments.find((seg) => seg.chromosome === d.chromosome);
        if (!segment) return null;

        const startAngle =
          segment.startAngle +
          (d.start / segment.length) * (segment.endAngle - segment.startAngle);
        const endAngle =
          segment.startAngle +
          (d.end / segment.length) * (segment.endAngle - segment.startAngle);

        const barArc = d3
          .arc()
          .innerRadius(innerRadius)
          .outerRadius(barScale(d.value))
          .startAngle(startAngle)
          .endAngle(endAngle);

        return barArc();
      })
      .attr("fill", (d) => {
        const segment = segments.find((seg) => seg.chromosome === d.chromosome);
        return segment ? segment.colour : "black";
      })
      .attr("stroke", "black");
  }, [canvasRef, config, data, segments]);

  return null;
};

export default Bar;
