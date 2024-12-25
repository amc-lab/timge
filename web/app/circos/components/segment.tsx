"use client";

import { useEffect } from "react";
import * as d3 from "d3";
import { Assembly, AssemblyConfig, GlobalConfig } from "../../types/genomes";

interface SegmentProps {
    data: {
        segments: Array<Assembly>;
        config: AssemblyConfig;
        globalConfig: GlobalConfig;
        divRef: any;
    };
    onSegmentsCreated?: (segmentsData: any[]) => void; // Callback to pass segment information to the parent
}

const Segment = ({data, onSegmentsCreated}: SegmentProps) => {
    const segments = data.segments;
    const config = data.config;
    const canvasRef = data.divRef;
    const globalConfig = data.globalConfig;

    const colorPalette = d3.scaleSequential(d3.interpolateSpectral);

    useEffect(() => {
        if (canvasRef.current && segments.length > 0) {
    
            let svg = d3.select(canvasRef.current).select("svg");
        
            if (svg.empty()) {
                svg = d3.select(canvasRef.current)
                    .append("svg")
                    .attr("width", globalConfig.canvasWidth)
                    .attr("height", globalConfig.canvasHeight);
            }

            const uniqueGroupClass = `group-1`;

            svg.selectAll(`g.${uniqueGroupClass}`).remove();

            const totalLength = segments.reduce(
                (acc, segment) => acc + (segment.end - segment.start),
                0
            );

            const chord = d3.chord()
                .padAngle(config.segmentPadding)
                .sortSubgroups(d3.descending);

            const arc = d3.arc()
                .innerRadius(config.segmentInnerRadius)
                .outerRadius(config.segmentOuterRadius);

            const axisGridlineArc = d3.arc()
                .innerRadius(config.segmentOuterRadius + config.segmentGridPadding)
                .outerRadius(config.segmentOuterRadius + config.segmentGridPadding);

            const matrix = segments.map((segment, i) =>
                segments.map((_, j) => (i === j ? segment.end - segment.start : 0))
            );

            const chords = chord(matrix);

            svg.selectAll("g").remove();

            const group = svg.append("g")
            .attr("class", uniqueGroupClass) // Add the unique class
            .attr("transform", `translate(${globalConfig.canvasWidth / 2}, ${globalConfig.canvasHeight / 2})`)
            .selectAll("g")
            .data(chords.groups)
            .join("g");

            group
            .append("path")
            .attr("fill", (d) => colorPalette(d.index / segments.length))
            .attr("d", arc)
            .attr("stroke", `${config.useStroke ? "black" : "none"}`)
            .on("mouseover", (event, d) => {
                d3.select(event.currentTarget).attr("filter", "brightness(0.95)");
            })
            .on("mouseout", (event, d) => {
                d3.select(event.currentTarget).attr("filter", "brightness(1)");
            })
            .append("title");

            const segmentData = chords.groups.map((d) => ({
                startAngle: d.startAngle,
                endAngle: d.endAngle,
                index: d.index,
                segment: segments[d.index],
                chromosome: segments[d.index].chromosome,
                length: segments[d.index].end - segments[d.index].start,
                colour: d3.color(colorPalette(d.index / segments.length)).formatHex(),
            }));

            if (onSegmentsCreated) {
                console.log(segmentData);
                onSegmentsCreated(segmentData);
            }

            group
                .append("text")
                .attr("transform", (d: { startAngle: number; endAngle: number; }) => {
                    const angle = (d.startAngle + d.endAngle) / 2;
                    const x = Math.sin(angle) * (config.segmentInnerRadius + (config.segmentOuterRadius - config.segmentInnerRadius) / 2);
                    const y = -Math.cos(angle) * (config.segmentInnerRadius + (config.segmentOuterRadius - config.segmentInnerRadius) / 2);
                    return `translate(${x}, ${y})`;
                })
                .attr("text-anchor", "middle")
                .attr("alignment-baseline", "middle")
                .attr("font-size", `10`)
                .attr("fill", "white")
                .text((d: { index: number; }) => segments[d.index].chromosome);

            const precision_vals = [100, 250, 500]
            const standardIntervalSizes = [500, 1000, 2500, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000]
            const expectedTickPrecision = totalLength / precision_vals[config.precision];
            const minorTickInterval = d3.bisect(standardIntervalSizes, expectedTickPrecision) === 0 ? standardIntervalSizes[0] : standardIntervalSizes[d3.bisect(standardIntervalSizes, expectedTickPrecision) - 1];
            const tickInterval = minorTickInterval * 5;

            if (config.showAxis) {
                group
                    .append("path")
                    .attr("d", axisGridlineArc)
                    .attr("stroke", "black")

                group
                    .append("g")
                    .selectAll("g")
                    .data((d: { index: string | number; startAngle: any; endAngle: any; }) => {
                    const segment = segments[d.index];
                    const scale = d3.scaleLinear()
                        .domain([segment.start, segment.end])
                        .range([d.startAngle, d.endAngle]);

                    const majorTicks = d3.range(segment.start, segment.end, tickInterval).map((value: number) => ({
                        value: value - segment.start == 0 ? null : d3.formatPrefix(",.0", 1e3)(value),
                        angle: scale(value),
                        isMajor: true,
                    }));

                    const minorTicks = d3.range(segment.start, segment.end, minorTickInterval).filter((value: number) => value % tickInterval !== 0).map((value: any) => ({
                        value: null,
                        angle: scale(value),
                        isMajor: false,
                    }));

                    return [...majorTicks, ...minorTicks];
                    })
                    .join("g")
                    .attr("transform", (d: { angle: number; }) => `rotate(${(d.angle * 180) / Math.PI - 90}) translate(${config.segmentOuterRadius + config.segmentGridPadding}, 0)`)
                    .call((g) =>
                    g
                        .append("line")
                        .attr("stroke", "currentColor")
                        .attr("x2", (d: { isMajor: any; }) => (d.isMajor ? config.tickLength * 2 : config.tickLength))
                    )
                    .call((g) =>
                    g
                        .filter((d: { isMajor: any; }) => d.isMajor)
                        .append("text")
                        .attr("x", 2 * config.tickLength + config.tickTextPadding)
                        .attr("dy", "0.35em")
                        .attr("font-size", `${config.axisLabelFontSize}`)
                        .attr("transform", (d: { angle: number; }) =>
                        d.angle > Math.PI ? `rotate(180) translate(-${4 * config.tickLength + 2 * config.tickTextPadding})` : null
                        )
                        .attr("text-anchor", (d: { angle: number; }) => (d.angle > Math.PI ? "end" : null))
                        .text((d: { value: { toLocaleString: () => any; }; }) => d.value?.toLocaleString())
                    );
            }
        }
    }, [segments, config, globalConfig]);

    return null;
};

export default Segment;
