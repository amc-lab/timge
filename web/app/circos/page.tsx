"use client";
import type { Assembly, Chord, AssemblyConfig, ChordConfig } from "../types/genomes";
import { useRef, useEffect, useState } from "react";
import * as d3 from "d3";
import FileUpload from "../../components/FileUpload";
import Segment from "./components/segment";
import Form from "../../components/Form";
import Chords from "./components/chord";

interface CircosProps {
    data: {
        segments: Array<Assembly>;
        chords: Array<Chord>;
    };
}

const defaultAssemblyConfig = {
    segmentPadding: 0.02,
    axisLabelFontSize: 10,
    showAxis: true,
    segmentInnerRadius: 200,
    segmentOuterRadius: 250,
    segmentGridPadding: 5,
    canvasWidth: 650,
    canvasHeight: 650,
    tickLength: 5,
    tickTextPadding: 2,
    precision: 1,
    useStroke: true,
};

const defaultChordConfig = {
    canvasWidth: 650,
    canvasHeight: 650,
};

const Circos = ({ data }: CircosProps) => {
    const canvasRef = useRef<HTMLDivElement>(null);
    const [segments, setSegments] = useState<Assembly[]>([]);
    const [chords, setChords] = useState<Chord[]>([]);
    const [assemblyConfig, setAssemblyConfig] = useState(defaultAssemblyConfig);
    const [chordConfig, setChordConfig] = useState(defaultChordConfig);

    const handleFileUpload = (file: File) => {
        const reader = new FileReader();
        reader.onload = (event) => {
            try {
                const json = JSON.parse(event.target?.result as string);
                setSegments(json);
            } catch (error) {
                console.error("Invalid JSON file:", error);
                alert("Uploaded file is not valid JSON.");
            }
        };
        reader.readAsText(file);
    };

    useEffect(() => {
        if (canvasRef.current) {
            const svg = d3.select(canvasRef.current).select("svg");

            if (svg.empty()) {
                d3.select(canvasRef.current)
                    .append("svg")
                    .attr("width", assemblyConfig.canvasWidth)
                    .attr("height", assemblyConfig.canvasHeight);
            }
        }
    }, [data, assemblyConfig.canvasWidth, assemblyConfig.canvasHeight]);

    const handleAssemblyConfigUpdate = (updatedConfig: Partial<AssemblyConfig>) => {
        setAssemblyConfig((prevConfig) => ({ ...prevConfig, ...updatedConfig }));
    };

    const handleChordConfigUpdate = (updatedConfig: Partial<ChordConfig>) => {
        setChordConfig((prevConfig) => ({ ...prevConfig, ...updatedConfig }));
    }

    return (
        <div>
            <FileUpload onFileUpload={handleFileUpload} />
            <div style={{ display: "flex", flexDirection: "row", height: "100vh" }}>
                <div style={{ flex: 5, display: "flex", alignItems: "center", justifyContent: "center" }}>
                    <div ref={canvasRef} style={{ border: "1px solid black" }}>
                        <Segment data={{ segments, config: assemblyConfig, divRef: canvasRef }} />
                        <Chords data={{ chords: chords, config: chordConfig, divRef: canvasRef }} />
                    </div>
                </div>

                <div style={{ flex: 5, display: "flex", alignItems: "center", justifyContent: "center" }}>
                    {segments.length > 0 && (
                        <Form handleAssemblyConfigUpdate={handleAssemblyConfigUpdate} defaultAssemblyConfig={assemblyConfig} 
                            handleChordConfigUpdate={handleChordConfigUpdate} defaultChordConfig={chordConfig} />
                    )}
                </div>
            </div>
        </div>
    );
};

export default Circos;
