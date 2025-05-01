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

                const arrowSpacing = 20;
                arcGroup.each(function(d, i) {
                    const currentInnerRadius = innerRadius + levels[i] * levelStep;
                    const currentOuterRadius = currentInnerRadius + trackWidth;
                    const midRadius = (currentInnerRadius + currentOuterRadius) / 2;
                
                    const startAngle = angleScale(d.chromStart);
                    const endAngle = angleScale(d.chromEnd);
                
                    const arcLength = midRadius * Math.abs(endAngle - startAngle);
                    const numArrows = Math.floor(arcLength / arrowSpacing);
                
                    const g = d3.select(this);
                
                    for (let j = 0; j < numArrows; j++) {
                        const t = (j + 0.5) / numArrows;
                        const angle = d.strand === "+" || d.strand === "."
                            ? startAngle + t * (endAngle - startAngle)
                            : endAngle - t * (endAngle - startAngle);
                
                        const x = midRadius * Math.cos(angle - Math.PI / 2);
                        const y = midRadius * Math.sin(angle - Math.PI / 2);
                
                        const arrowSize = 6;
                
                        let angleDeg = angle * 180 / Math.PI + 90;
                
                        if (d.strand === "-") {
                            angleDeg += 180;
                        }
                
                        g.append("polygon")
                            .attr("points", `${-arrowSize / 2},${arrowSize} ${arrowSize / 2},${arrowSize} 0,0`)
                            .attr("fill", "white")
                            .attr("transform", `translate(${x},${y}) rotate(${angleDeg})`);
                    }
                });

                arcGroup.append("text")
                .text(d => d.name)
                .attr("font-family", "Arial")
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
                })
                .attr("display", (d, i) => {
                    // if annotation is too close to the next one, hide it
                    const nextAnnotation = sortedAnnotations[i + 1];
                    if (nextAnnotation) {
                        const nextStartAngle = angleScale(nextAnnotation.chromStart);
                        const nextEndAngle = angleScale(nextAnnotation.chromEnd);
                        const currentStartAngle = angleScale(d.chromStart);
                        const currentEndAngle = angleScale(d.chromEnd);

                        const nextMidAngle = (nextStartAngle + nextEndAngle) / 2;
                        const currentMidAngle = (currentStartAngle + currentEndAngle) / 2;

                        const angleDiff = Math.abs(nextMidAngle - currentMidAngle);
                        return angleDiff < Math.PI / 20 ? "none" : "visible";
                    }
                    return "block";
                });
            
        });

    }, [annotations, config, segments, globalConfig]);

    return null;
};

export default Annotation;
