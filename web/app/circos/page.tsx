"use client";
import { type Assembly, type Chord, type AssemblyConfig, type ChordConfig, type GlobalConfig } from "../types/genomes";
import { useRef, useEffect, useState } from "react";
import * as d3 from "d3";
import FileUpload from "../../components/FileUpload";
import Form from "../../components/Form";
import Tracks from "./tracks";
import { Track, TrackType } from "./config/track";
import { sampleChords } from "./config/sampleChords";
import { defaultAssemblyConfig, defaultChordConfig, defaultGlobalConfig } from "./config/defaultConfigs";

interface CircosProps {
    data: {
        segments: Array<Assembly>;
        chords: Array<Chord>;
    };
}

const Circos = ({ data }: CircosProps) => {
    const canvasRef = useRef<HTMLDivElement>(null);
    const [segments, setSegments] = useState<Assembly[]>([]);
    const [chords, setChords] = useState<Chord[]>(sampleChords);
    const [assemblyConfig, setAssemblyConfig] = useState(defaultAssemblyConfig);
    const [chordConfig, setChordConfig] = useState(defaultChordConfig);
    const [globalConfig, setGlobalConfig] = useState(defaultGlobalConfig);
    const [tracks, setTracks] = useState<Track[]>([]);

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
                    .attr("width", globalConfig.canvasWidth)
                    .attr("height", globalConfig.canvasHeight)
            }
        }
    }, [data, globalConfig]);

    const handleAssemblyConfigUpdate = (updatedConfig: Partial<AssemblyConfig>) => {
        setAssemblyConfig((prevConfig) => ({ ...prevConfig, ...updatedConfig }));
        setChordConfig((prevConfig) => ({ ...prevConfig, outerRadius: updatedConfig.segmentInnerRadius ?? prevConfig.outerRadius }));
    };

    const handleGlobalConfigUpdate = (updatedConfig: Partial<GlobalConfig>) => {
        setGlobalConfig((prevConfig) => ({ ...prevConfig, ...updatedConfig }));
    }

    const handleTrackConfigUpdate = (index: number, updatedConfig: any) => {
        setTracks((prevTracks) => {
            const newTracks = [...prevTracks];
            newTracks[index] = { ...newTracks[index], config: updatedConfig };
            return newTracks;
        });
    }

    useEffect(() => {
        setTracks([
            {
                trackType: TrackType.Karyotype,
                data: { segments, globalConfig, divRef: canvasRef },
                config: assemblyConfig,
            },
            {
                trackType: TrackType.Chord,
                data: { chords, globalConfig, divRef: canvasRef },
                config: chordConfig,
            },
        ]);
    }, [segments, chords, assemblyConfig, chordConfig, globalConfig]);

    return (
        <div>
            <FileUpload onFileUpload={handleFileUpload} />
            <div style={{ display: "flex", flexDirection: "row", height: "100vh" }}>
                <div style={{ flex: 5, display: "flex", alignItems: "center", justifyContent: "center" }}>
                    <div ref={canvasRef} style={{ border: "1px solid black" }}>
                        <Tracks tracks={tracks} />
                    </div>
                </div>

                <div style={{ flex: 5, display: "flex", alignItems: "center", justifyContent: "center" }}>
                    {segments.length > 0 && (
                        <Form 
                            tracks = {tracks}
                            handleTrackConfigUpdate={handleTrackConfigUpdate}
                            handleAssemblyConfigUpdate={handleAssemblyConfigUpdate} defaultAssemblyConfig={assemblyConfig} 
                            defaultChordConfig={chordConfig}
                            handleGlobalConfigUpdate={handleGlobalConfigUpdate} defaultGlobalConfig={globalConfig}
                        />
                    )}
                </div>
            </div>
        </div>
    );
};

export default Circos;
