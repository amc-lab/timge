"use client";
import type { Assembly, Chord, AssemblyConfig, ChordConfig } from "../types/genomes";
import { useRef, useEffect, useState } from "react";
import * as d3 from "d3";
import FileUpload from "../../components/FileUpload";
import Segment from "./segment";
import Form from "../../components/Form";

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

const defaultChordConfig = {};

const Circos = ({ data }: CircosProps) => {
    const canvasRef = useRef<HTMLDivElement>(null);
    const [segments, setSegments] = useState<Assembly[]>([]);
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
            const svg = d3.select(canvasRef.current)
                .append("svg")
                .attr("width", assemblyConfig.canvasWidth)
                .attr("height", assemblyConfig.canvasHeight);

            svg.call(d3.create);
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
                {/* Left Section */}
                <div style={{ flex: 5, display: "flex", alignItems: "center", justifyContent: "center" }}>
                    <Segment data={{ segments, config: assemblyConfig }} />
                </div>

                {/* Right Section (Form) */}
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
