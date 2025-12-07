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
  segments: Array<any>;
  idx: number;
}

const Ring = ({ data, config, segments, idx }: RingProps) => {
  const { sequences, divRef, globalConfig } = data;

  useEffect(() => {
    if (!divRef.current || sequences.length === 0 || segments.length === 0 || config.hide ) return;

    let svg = d3.select(divRef.current).select("svg");
    if (svg.empty()) {
      d3.select(divRef.current).select("canvas").remove(); // Remove old canvas if exists
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
      console.log(segments);
      // const angle = (segmentLength / totalLength) * 2 * Math.PI;
      const angle = (segmentLength / segments[i].length) * (segments[i].endAngle - segments[i].startAngle);

      const arc = d3
        .arc()
        .innerRadius(config.innerRadius)
        .outerRadius(config.innerRadius + config.trackWidth)
        .startAngle(startAngle)
        .endAngle(startAngle + angle);

      const axisGridlineArc = d3
        .arc()
        .innerRadius(
          config.innerRadius +
          config.trackWidth +
          config.gridPadding,
        )
        .outerRadius(
          config.innerRadius +
          config.trackWidth +
          config.gridPadding,
        )
        .startAngle(startAngle)
        .endAngle(startAngle + angle);

      group
        .append("path")
        .attr("d", arc as any)
        .attr("fill", d3.interpolateSpectral(idx / 5))
        .attr("stroke", "black")
        .on("mouseover", (event, d) => {
          d3.select(event.currentTarget).attr("filter", "brightness(0.95)");
        })
        .on("mouseout", (event, d) => {
          d3.select(event.currentTarget).attr("filter", "brightness(1)");
        });

        group
        .append("text")
        .attr("transform", (d: { startAngle: number; endAngle: number }) => {
            const midAngle = startAngle + angle / 2;
            const x = Math.sin(midAngle) * (config.innerRadius + config.trackWidth / 2);
            const y = -Math.cos(midAngle) * (config.innerRadius + config.trackWidth / 2);
            return `translate(${x}, ${y})`;
        })
        .attr("text-anchor", "middle")
        .attr("alignment-baseline", "middle")
        .attr("font-family", "Arial")
        .attr("font-size", 10)
        .attr("fill", "black")
        .text((d: { index: number }) => segment.chromosome);

      const scale = d3
        .scaleLinear()
        .domain([segment.start, segment.end])
        .range([startAngle, startAngle + angle]);


      if (config.showAxis) {
        group.append("path").attr("d", axisGridlineArc).attr("stroke", "black");

        const precision_vals = [100, 250, 500];
        const standardIntervalSizes = [
          50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000, 100000, 250000, 500000,
          1000000,
        ];
        const expectedTickPrecision =
            totalLength / precision_vals[config.precision];
        const minorTickInterval =
            d3.bisect(standardIntervalSizes, expectedTickPrecision) === 0
            ? standardIntervalSizes[0]
            : standardIntervalSizes[
                d3.bisect(standardIntervalSizes, expectedTickPrecision) - 1
                ];
        const tickInterval = minorTickInterval * 5;

        const ticks = d3.range(segment.start, segment.end, tickInterval);
        const tickGroup = group.append("g");

        tickGroup
            .selectAll("g")
            .data(() => {
            const majorTicks = d3
                .range(segment.start, segment.end, tickInterval)
                .map((value) => {
                const prefixScale =
                    config.metricPrefix === "k"
                    ? 1e3
                    : config.metricPrefix === "M"
                    ? 1e6
                    : config.metricPrefix === "G"
                    ? 1e9
                    : 1;

                const formattedValue =
                    value - segment.start === 0
                    ? null
                    : d3.formatPrefix(".2", prefixScale)(value);

                return {
                    value: formattedValue,
                    angle: scale(value),
                    isMajor: true,
                };
                });

            const minorTicks = d3
                .range(segment.start, segment.end, minorTickInterval)
                .filter((value: number) => value % tickInterval !== 0)
                .map((value: any) => ({
                value: null,
                angle: scale(value),
                isMajor: false,
                }));

            return [...majorTicks, ...minorTicks];
            })
            .join("g")
            .attr(
            "transform",
            (d: { angle: number }) =>
                `rotate(${(d.angle * 180) / Math.PI - 90}) translate(${config.innerRadius + config.trackWidth + config.gridPadding}, 0)`
            )
            .call((g) =>
            g
                .append("line")
                .attr("stroke", "currentColor")
                .attr("x2", (d: { isMajor: any }) =>
                d.isMajor ? config.tickLength * 2 : config.tickLength
                )
            )
            .call((g) =>
            g
                .filter((d: { isMajor: any }) => d.isMajor)
                .append("text")
                .attr("x", 2 * config.tickLength + config.tickTextPadding)
                .attr("dy", "0.35em")
                .attr("font-size", `${config.axisLabelFontSize}`)
                .attr("transform", (d: { angle: number }) =>
                d.angle > Math.PI
                    ? `rotate(180) translate(-${4 * config.tickLength + 2 * config.tickTextPadding})`
                    : null
                )
                .attr("text-anchor", (d: { angle: number }) =>
                d.angle > Math.PI ? "end" : null
                )
                .text((d: { value: { toLocaleString: () => any } }) =>
                d.value?.toLocaleString()
                )
            );
          }

      startAngle += angle;
    });
  }, [sequences, config, segments, globalConfig]);

  return null;
};

export default Ring;
