import { useEffect } from "react";
import * as d3 from "d3";
import { GlobalConfig } from "@/app/types/genomes";

interface HighlightProps {
    segmentStartIdx: number;
    segmentEndIdx: number;
    segmentStartPos: number;
    segmentEndPos: number;
    divRef: any;
    segments: Array<any>;
    globalConfig: GlobalConfig;
    innerRadius: number;
    width: number;
}

const Highlight = ({ globalConfig, segments, segmentStartIdx, segmentEndIdx, segmentStartPos, segmentEndPos, divRef, innerRadius, width }: HighlightProps) => {
    console.log(segmentStartIdx, segmentEndIdx, segmentStartPos, segmentEndPos, divRef);
    useEffect(() => {
        console.log("Highlighting", segments);
        if (!divRef || !divRef.current || segments.length === 0) return;

        let svg = d3.select(divRef.current).select("svg");
    
        const uniqueGroupClass = `highlight-group`;
        svg.selectAll(`g.${uniqueGroupClass}`).remove();
    
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
        
    }, [segmentStartIdx, segmentEndIdx, segmentStartPos, segmentEndPos, divRef, segments]);
    
    return null;
};

export default Highlight;