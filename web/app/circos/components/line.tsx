import { GlobalConfig, LineConfig, LineData } from "@/app/types/genomes";
import { useEffect } from "react";
import * as d3 from "d3";

interface LineProps {
    data: {
        values: Array<LineData>;
        globalConfig: GlobalConfig;
        divRef: any;
    };
    config: LineConfig;
    segments: Array<any>;
    idx: number;
}

const Line = ({ data, config, segments, idx }: LineProps) => {
    const { values, divRef, globalConfig } = data;

    console.log("Line component props", data, config, segments, idx);

    useEffect(() => {
        if (!divRef.current || values.length === 0 || segments.length === 0 || config.hide) return;
    
        let svg = d3.select(divRef.current).select("svg");
        if (svg.empty()) {
            d3.select(divRef.current).select("canvas").remove();
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
            .attr("transform", `translate(${globalConfig.canvasWidth / 2}, ${globalConfig.canvasHeight / 2})`);
    
        segments.forEach((segment) => {
            console.log(segment, segment.chrom);
            const segmentValues = values.filter(d => d.chrom === segment.chromosome);
            console.log(segmentValues);
            const max = d3.max(segmentValues, d => d.value);
            const min = d3.min(segmentValues, d => d.value);
            console.log(segmentValues);
            console.log(max);
            console.log(min);
            segmentValues.forEach((d) => {
                d.value = d.value / max * config.trackWidth;
            });
            console.log(segmentValues);

            const angleScale = d3.scaleLinear()
                .domain([0, segment.length])
                .range([segment.startAngle, segment.endAngle]);

            const radialLine = d3.lineRadial()
                .angle(d => angleScale(d.chromStart))
                .radius(d => config.innerRadius + d.value)
                .curve(d3.curveBasis);

            group.append("path")
                .datum(segmentValues)
                .attr("d", radialLine)
                .attr("fill", "none")
                .attr("stroke", config.colour || "black")
                .attr("stroke-width", 1);
        });
    }, [values, config, segments, globalConfig]);

    return null;
};

export default Line;
