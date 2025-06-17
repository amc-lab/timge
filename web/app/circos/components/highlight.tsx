import { useEffect } from "react";
import * as d3 from "d3";
import { GlobalConfig, HighlightConfig } from "@/app/types/genomes";

interface HighlightProps {
    data: {
        // globalConfig: GlobalConfig;
        divRef: any;
    }
    config: HighlightConfig;
    segments: Array<any>;
    dependencies?: any;
    globalConfig?: GlobalConfig;
}

const Highlight = ({  data, config, segments, dependencies, globalConfig }: HighlightProps) => {
    const { divRef } = data;

    useEffect(() => {
        if (!divRef || !divRef.current || segments.length === 0) return;

        let svg = d3.select(divRef.current).select("svg");
    
        const uniqueGroupClass = `highlight-group`;
        svg.selectAll(`g.${uniqueGroupClass}`).remove();

        if (!globalConfig.showHighlight)
            return;

        const { innerRadius, width } = config;
        if (!dependencies) {
            dependencies = { chr: undefined, start: undefined, end: undefined };
        }

        const { chr, start, end } = dependencies;
        if (!chr || chr === "all" || start === undefined || end === undefined) return;

        const segmentStartIdx = segments.findIndex(segment => segment.chromosome === chr);
        const segmentEndIdx = segments.findIndex(segment => segment.chromosome === chr);
        const segmentStartPos = start;
        const segmentEndPos = end;      
    
        const group = svg
            .append("g")
            .attr("class", uniqueGroupClass)
            .attr("transform", `translate(${globalConfig.canvasWidth / 2}, ${globalConfig.canvasHeight / 2})`);
    
        const startAngle = segments[segmentStartIdx].startAngle + 
            segmentStartPos / segments[segmentStartIdx].length * (segments[segmentStartIdx].endAngle - segments[segmentStartIdx].startAngle);
        const endAngle = segments[segmentEndIdx].startAngle +
            segmentEndPos / segments[segmentEndIdx].length * (segments[segmentEndIdx].endAngle - segments[segmentEndIdx].startAngle);

        const highlightArc = d3.arc()
            .innerRadius(innerRadius - 5)
            .outerRadius(innerRadius + width + 5)
            .startAngle(startAngle)
            .endAngle(endAngle);

        group.append("path")
            .attr("d", highlightArc as any)
            .attr("fill", "blue")
            .attr("opacity", 0.1);
        
    }, [config, divRef, segments, dependencies, globalConfig]);
    
    return null;
};

export default Highlight;