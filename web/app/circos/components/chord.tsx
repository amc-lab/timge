import React, { useEffect } from 'react';
import { Chord, ChordConfig, GlobalConfig } from '@/app/types/genomes';
import { useRef } from 'react';
import * as d3 from 'd3';

interface ChordProps {
    data: {
        chords: Array<Chord>;
        config: ChordConfig;
        globalConfig: GlobalConfig;
        divRef: any;
    };
}

const Chords = ({ data }: ChordProps) => {
    const canvasRef = data.divRef;
    const chords = data.chords;
    const config = data.config;
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

        const group = svg.append("g")
            .attr("transform", `translate(${globalConfig.canvasWidth / 2}, ${globalConfig.canvasHeight / 2})`);

        const radius = 200;

        group.selectAll("path")
            .data(chords)
            .join("path")
            .attr("d", d3.ribbon()
                .radius(radius)
                .source(d => ({ startAngle: (d.source_start / radius) * 2 * Math.PI, endAngle: (d.source_end / radius) * 2 * Math.PI }))
                .target(d => ({ startAngle: (d.target_start / radius) * 2 * Math.PI, endAngle: (d.target_end / radius) * 2 * Math.PI }))
            )
            .attr("fill", "none")
            .attr("stroke", "#000")
            .attr("stroke-width", 1);

        console.log(chords);
    }, [canvasRef, chords, config, globalConfig]);

    return null;
};

export default Chords;
