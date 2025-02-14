import React, { useEffect } from "react";
import * as d3 from "d3";
import { GlobalConfig, RingConfig, RingData } from "@/app/types/genomes";

interface RingProps {
  data: {
    sequences: Array<RingData>;
    globalConfig: GlobalConfig;
    divRef: any;
  };
  config: RingConfig;
  idx: number;
}

const Ring = ({ data, config, idx }: RingProps) => {
  const { sequences, divRef, globalConfig } = data;

  useEffect(() => {
    if (!divRef.current || sequences.length === 0) return;

    let svg = d3.select(divRef.current).select("svg");
    if (svg.empty()) {
      svg = d3
        .select(divRef.current)
        .append("svg")
        .attr("width", globalConfig.canvasWidth)
        .attr("height", globalConfig.canvasHeight);
    }

    const uniqueGroupClass = `group-${idx}`;
    svg.selectAll(`g.${uniqueGroupClass}`).remove();

    const group = svg
      .append("g")
      .attr("class", uniqueGroupClass)
      .attr(
        "transform",
        `translate(${globalConfig.canvasWidth / 2}, ${globalConfig.canvasHeight / 2})`
      );

    const totalLength = sequences.reduce(
      (acc, segment) => acc + (segment.end - segment.start),
      0
    );

    let startAngle = 0;
    sequences.forEach((segment, i) => {
      const segmentLength = segment.end - segment.start;
      const angle = (segmentLength / totalLength) * 2 * Math.PI;

      const arc = d3
        .arc()
        .innerRadius(config.innerRadius)
        .outerRadius(config.innerRadius + config.trackWidth)
        .startAngle(startAngle)
        .endAngle(startAngle + angle);

      group
        .append("path")
        .attr("d", arc as any)
        .attr("fill", d3.interpolateSpectral(i / sequences.length))
        .attr("stroke", "black");

      const axisGroup = group.append("g").attr("class", `axis-${i}`);

      startAngle += angle;
    });
  }, [sequences, config, globalConfig]);

  return null;
};

export default Ring;
