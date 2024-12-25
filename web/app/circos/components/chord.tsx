import React, { useEffect } from 'react';
import { Chord, ChordConfig } from '@/app/types/genomes';
import { useRef } from 'react';
import * as d3 from 'd3';

interface ChordProps {
    data: {
        chords: Array<Chord>;
        config: ChordConfig;
        divRef: any;
    };
}

const Chords = ({data}: ChordProps) => {
    const canvasRef = data.divRef;
    const chords = data.chords;
    const config = data.config;

    useEffect(() => {
        if (!canvasRef.current) return;
    
        let svg = d3.select(canvasRef.current).select("svg");
    
        if (svg.empty()) {
            svg = d3.select(canvasRef.current)
                .append("svg")
                .attr("width", config.canvasWidth)
                .attr("height", config.canvasHeight);
        }
    
        // svg.append("circle")
        //     .attr("cx", config.canvasWidth / 2)
        //     .attr("cy", config.canvasHeight / 2)
        //     .attr("r", 200)
        //     .attr("fill", "red");
        
    }, [canvasRef, chords, config]);
    

    return null;
}

export default Chords;