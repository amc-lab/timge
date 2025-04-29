"use client";
import { AnnotationConfig, AnnotationData, GlobalConfig } from "@/app/types/genomes";
import { useEffect } from "react";
import * as d3 from "d3";

interface AnnotationProps {
    data: {
        annotations: Array<AnnotationData>;
        globalConfig: GlobalConfig;
        divRef: any;
    };
    config: AnnotationConfig;
    segments: Array<any>;
    idx: number;
}

const Annotation = ({ data, config, segments, idx }: AnnotationProps) => {
    const { annotations, divRef, globalConfig } = data;

    useEffect(() => {
        if (!divRef.current || annotations.length === 0 || segments.length === 0 || config.hide) return;

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

        const trackWidth = config.trackWidth;
        const trackPadding = config.trackPadding || 5;
        const textPadding = config.textPadding || 5;
        const textFontSize = config.textFontSize || 10;
        const innerRadius = config.innerRadius;

        const levelStep = trackWidth + trackPadding + textFontSize + textPadding;

        segments.forEach((segment) => {
            const _annotations = annotations.filter(d => d.chrom === segment.chromosome);

            const angleScale = d3.scaleLinear()
                .domain([0, segment.length])
                .range([segment.startAngle, segment.endAngle]);

            const sortedAnnotations = [..._annotations].sort((a, b) => a.chromStart - b.chromStart);

            const levels: number[] = [];
            sortedAnnotations.forEach((d, i) => {
                let level = 0;
                const dStart = angleScale(d.chromStart);
                const dEnd = angleScale(d.chromEnd);

                while (true) {
                    const conflict = sortedAnnotations.some((other, j) => {
                        if (i === j) return false;
                        if (levels[j] !== level) return false;
                        const oStart = angleScale(other.chromStart);
                        const oEnd = angleScale(other.chromEnd);
                        return !(dEnd <= oStart || dStart >= oEnd);
                    });
                    if (!conflict) break;
                    level += 1;
                }
                levels[i] = level;
            });

            const arcGroup = group.selectAll(`.annotation-arc-${segment.chromosome}`)
                .data(sortedAnnotations)
                .enter()
                .append("g")
                .attr("class", `annotation-arc-${segment.chromosome}`);

            arcGroup.append("path")
                .attr("d", (d, i) => {
                    const currentInnerRadius = innerRadius + levels[i] * levelStep;
                    const currentOuterRadius = currentInnerRadius + trackWidth;

                    const arc = d3.arc<AnnotationData>()
                        .innerRadius(currentInnerRadius)
                        .outerRadius(currentOuterRadius)
                        .startAngle(() => angleScale(d.chromStart))
                        .endAngle(() => angleScale(d.chromEnd));

                    return arc(d);
                })
                .attr("fill", d => d.itemRgb ? `rgb(${d.itemRgb})` : "steelblue")
                .attr("stroke", "white")
                .attr("stroke-width", 0.5);

                arcGroup.append("text")
                .text(d => d.name)
                .attr("font-size", textFontSize)
                .attr("fill", d => d.itemRgb ? `rgb(${d.itemRgb})` : "black")
                .attr("text-anchor", "middle")
                .attr("alignment-baseline", "middle")
                .attr("transform", (d, i) => {
                    const midAngle = (angleScale(d.chromStart) + angleScale(d.chromEnd)) / 2;
            
                    const currentInnerRadius = innerRadius + levels[i] * levelStep;
                    const labelRadius = currentInnerRadius - textPadding - (textFontSize / 2);
            
                    const x = labelRadius * Math.cos(midAngle - Math.PI / 2);
                    const y = labelRadius * Math.sin(midAngle - Math.PI / 2);
                    return `translate(${x},${y}) rotate(${(midAngle * 180 / Math.PI)})`;
                });
            
        });

    }, [annotations, config, segments, globalConfig]);

    return null;
};

export default Annotation;
