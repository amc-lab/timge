import { useEffect } from 'react';
import { Chord, ChordConfig, GlobalConfig } from '@/app/types/genomes';
import * as d3 from 'd3';

interface ChordProps {
    data: {
        chords: Array<Chord>;
        globalConfig: GlobalConfig;
        divRef: any;
    };
    config: ChordConfig;
    segments: Array<any>;
}

const Chords = ({ data, config, segments }: ChordProps) => {
    const canvasRef = data.divRef;
    const chords = data.chords;
    const globalConfig = data.globalConfig;

    useEffect(() => {
        if (!canvasRef.current) return;
        let svg = d3.select(canvasRef.current).select("svg");

        if (svg.empty()) {
            svg = d3.select(canvasRef.current)
                .append("svg")
                .attr("width", globalConfig.canvasWidth)
                .attr("height", globalConfig.canvasHeight);
        }

        const uniqueGroupClass = `group-2`;

        svg.selectAll(`g.${uniqueGroupClass}`).remove();

        const group = svg.append("g")
            .attr("class", uniqueGroupClass)
            .attr("transform", `translate(${globalConfig.canvasWidth / 2}, ${globalConfig.canvasHeight / 2})`);

        const radius = config.outerRadius;
        const chord_padding = config.chordPadding;
        const chord_radius = radius - chord_padding;

        group.selectAll("path")
            .data(chords)
            .join("path")
            .attr("d", d3.ribbon()
                .radius(chord_radius)
                .source(d => (
                    { 
                        startAngle: segments.find(segment => segment.chromosome === d.source_chromosome)?.startAngle
                        + (d.source_start / (segments.find(segment => segment.chromosome === d.source_chromosome)?.length)) * (segments.find(segment => segment.chromosome === d.source_chromosome)?.endAngle - segments.find(segment => segment.chromosome === d.source_chromosome)?.startAngle), 
                        endAngle: segments.find(segment => segment.chromosome === d.source_chromosome)?.startAngle
                        + (d.source_end / (segments.find(segment => segment.chromosome === d.source_chromosome)?.length)) * (segments.find(segment => segment.chromosome === d.source_chromosome)?.endAngle - segments.find(segment => segment.chromosome === d.source_chromosome)?.startAngle),
                    }
                ))
                .target(d => ({ 
                        startAngle: segments.find(segment => segment.chromosome === d.target_chromosome)?.startAngle
                        + (d.target_start / (segments.find(segment => segment.chromosome === d.target_chromosome)?.length)) * (segments.find(segment => segment.chromosome === d.target_chromosome)?.endAngle - segments.find(segment => segment.chromosome === d.target_chromosome)?.startAngle), 
                        endAngle: segments.find(segment => segment.chromosome === d.target_chromosome)?.startAngle
                        + (d.target_end / (segments.find(segment => segment.chromosome === d.target_chromosome)?.length)) * (segments.find(segment => segment.chromosome === d.target_chromosome)?.endAngle - segments.find(segment => segment.chromosome === d.target_chromosome)?.startAngle), 
                 }))
            )
            .attr("fill", d => {
                const sourceSegment = segments.find(segment => segment.chromosome === d.source_chromosome);
                return sourceSegment?.colour || "gray";
            })
            .attr("opacity", config.opacity)
            .attr("stroke", `${config.useStroke ? "black" : "none"}`)
            .on("mouseover", function () {
                const currentOpacity = parseFloat(d3.select(this).attr("opacity"));
                d3.select(this).attr("opacity", Math.min(currentOpacity + 0.3, 1)); // Increase opacity but cap at 1
            })
            .on("mouseout", function () {
                const currentOpacity = parseFloat(d3.select(this).attr("opacity"));
                d3.select(this).attr("opacity", Math.max(currentOpacity - 0.3, 0)); // Decrease opacity but ensure it doesn't go below 0
            });

    }, [canvasRef, chords, config, globalConfig]);

    return null;
};

export default Chords;
