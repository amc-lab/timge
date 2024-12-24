"use client";
import type { Assembly, Chord, AssemblyConfig } from "../types/genomes";
import { useRef, useEffect, useState } from "react";
import * as d3 from "d3";
import FileUpload from "../../components/FileUpload";
import Segment from "./segment";
import {AssemblyForm} from "../../components/Form";

interface CircosProps {
    data: {
        segments: Array<Assembly>;
        chords: Array<Chord>;
    };
}

const defaultConfig = {
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

const Circos = ({ data }: CircosProps) => {
    const canvasRef = useRef<HTMLDivElement>(null);
    const [segments, setSegments] = useState<Assembly[]>([]);
    const [config, setConfig] = useState(defaultConfig);

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
                .attr("width", config.canvasWidth)
                .attr("height", config.canvasHeight);

            svg.call(d3.create);
        }
    }, [data, config.canvasWidth, config.canvasHeight]);

    const handleConfigUpdate = (updatedConfig: Partial<AssemblyConfig>) => {
        setConfig((prevConfig) => ({ ...prevConfig, ...updatedConfig }));
    };

    return (
        <>
            <FileUpload onFileUpload={handleFileUpload} />
            <Segment data={{ segments, config }} />
            {segments.length > 0 && (
                <AssemblyForm
                    onUpdate={handleConfigUpdate}
                    defaultConfig={defaultConfig}
                />
            )}
        </>
    );
};

export default Circos;
