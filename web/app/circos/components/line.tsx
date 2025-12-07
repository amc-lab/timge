import { GlobalConfig, LineConfig, LineData } from "@/app/types/genomes";
import { useEffect } from "react";
import * as d3 from "d3";
import { View } from "@/store/features/views/types";

interface LineProps {
    data: {
        values: Array<LineData>;
        globalConfig: GlobalConfig;
        divRef: any;
    };
    config: LineConfig;
    segments: Array<any>;
    idx: number;
    trackName?: string;
}

const Line = ({ data, config, segments, idx, trackName }: LineProps) => {
    const { values, divRef, globalConfig } = data;

    console.log("Line component props", data, config, segments, idx);

    useEffect(() => {
        if (!divRef.current || values.length === 0 || segments.length === 0 || config.hide) return;
    
        let svg = d3.select(divRef.current).select("svg");
        // if (svg.empty()) {
        //     d3.select(divRef.current).select("canvas").remove();
        //     svg = d3
        //         .select(divRef.current)
        //         .append("svg")
        //         .attr("width", globalConfig.canvasWidth)
        //         .attr("height", globalConfig.canvasHeight);
        // }

                if (svg.empty()) {
                  d3.select(divRef.current)
                    .append("svg")
                    .attr("viewBox", `0 0 ${globalConfig.canvasWidth} ${globalConfig.canvasHeight}`)
                    .attr("preserveAspectRatio", "xMidYMid meet")
                    .style("width", "100%")      // fill the container’s width
                    .style("height", "auto")     // height scaled to preserve aspect
                    .style("display", "block");  // remove inline-gap artifacts
                }

        const uniqueGroupClass = `group-${idx}`;
        svg.selectAll(`g.${uniqueGroupClass}`).remove();
    
        const group = svg
            .append("g")
            .attr("class", uniqueGroupClass)
            .attr("transform", `translate(${globalConfig.canvasWidth / 2}, ${globalConfig.canvasHeight / 2})`);
    
        segments.forEach((segment) => {
            const segLen = segment.length;
            const byPos = new Map<number, LineData>(values
                            .filter(d => d.chrom === segment.chromosome)
                            .map(d => [d.chromStart, d])
                    );
            const fullValues: Array<LineData | undefined> = Array.from({ length: segLen }, (_, i) =>
            byPos.get(i)
            );

            const validValues = fullValues.filter(d => d);
            const max = d3.max(validValues, d => d!.value)!;
            validValues.forEach(d => { if (d) d.value = (d.value / max) * config.trackWidth; });

            const angleScale = d3.scaleLinear()
            .domain([0, segLen])
            .range([segment.startAngle, segment.endAngle]);

            const radialLine = d3.lineRadial<LineData>()
            .defined(d => d !== undefined)
            .angle(d => angleScale(d!.chromStart))
            .radius(d => config.innerRadius + d!.value)
            .curve(d3.curveBasis);

            const runs: LineData[][] = [];
            let run: LineData[] = [];
            fullValues.forEach((d) => {
            if (d) run.push(d);
            else if (run.length) { runs.push(run); run = []; }
            });
            if (run.length) runs.push(run);

            runs.forEach(run => {
            group.append("path")
                .datum(run)
                .attr("d", radialLine)
                .attr("fill", "none")
                .attr("stroke", config.colour || "black")
                .attr("stroke-width", 1);
            });
        });
    }, [values, config, segments, globalConfig]);

    return null;
};

export default Line;
