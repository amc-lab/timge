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
    tickLength: 5,
    tickTextPadding: 2,
    precision: 1,
    useStroke: true,
};

const defaultChordConfig = {
};

const defaultGlobalConfig = {
    canvasWidth: 650,
    canvasHeight: 650,
};

const Circos = ({ data }: CircosProps) => {
    const canvasRef = useRef<HTMLDivElement>(null);
    const [segments, setSegments] = useState<Assembly[]>([]);
    const [chords, setChords] = useState<Chord[]>([]);
    const [assemblyConfig, setAssemblyConfig] = useState(defaultAssemblyConfig);
    const [chordConfig, setChordConfig] = useState(defaultChordConfig);
    const [globalConfig, setGlobalConfig] = useState(defaultGlobalConfig);

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

    const chordData = [
        {
            "source_chr": "os10",
            "source_start": 727761,
            "source_end": 862195,
            "target_chr": "os13",
            "target_start": 503711,
            "target_end": 727521
        },
        {
            "source_chr": "os11",
            "source_start": 53706,
            "source_end": 119898,
            "target_chr": "os12",
            "target_start": 132510,
            "target_end": 318392
        },
        {
            "source_chr": "os14",
            "source_start": 890047,
            "source_end": 1176200,
            "target_chr": "os9",
            "target_start": 709549,
            "target_end": 839923
        },
        {
            "source_chr": "os3",
            "source_start": 533225,
            "source_end": 603174,
            "target_chr": "os13",
            "target_start": 589625,
            "target_end": 824058
        },
        {
            "source_chr": "os8",
            "source_start": 623161,
            "source_end": 674225,
            "target_chr": "os4",
            "target_start": 237627,
            "target_end": 344041
        },
        {
            "source_chr": "os1",
            "source_start": 171532,
            "source_end": 213241,
            "target_chr": "os4",
            "target_start": 557498,
            "target_end": 632494
        },
        {
            "source_chr": "os7",
            "source_start": 611199,
            "source_end": 688226,
            "target_chr": "os14",
            "target_start": 1369178,
            "target_end": 1386275
        },
        {
            "source_chr": "os5",
            "source_start": 615607,
            "source_end": 711831,
            "target_chr": "os3",
            "target_start": 326396,
            "target_end": 406051
        },
        {
            "source_chr": "os6",
            "source_start": 337855,
            "source_end": 429067,
            "target_chr": "os8",
            "target_start": 241697,
            "target_end": 360185
        },
        {
            "source_chr": "os3",
            "source_start": 508749,
            "source_end": 511939,
            "target_chr": "os5",
            "target_start": 146236,
            "target_end": 179771
        },
        {
            "source_chr": "os4",
            "source_start": 313430,
            "source_end": 414960,
            "target_chr": "os14",
            "target_start": 1275583,
            "target_end": 1347285
        },
        {
            "source_chr": "os7",
            "source_start": 518785,
            "source_end": 636924,
            "target_chr": "os1",
            "target_start": 180067,
            "target_end": 180884
        },
        {
            "source_chr": "os4",
            "source_start": 56690,
            "source_end": 83076,
            "target_chr": "os14",
            "target_start": 234382,
            "target_end": 325789
        },
        {
            "source_chr": "os4",
            "source_start": 515704,
            "source_end": 593275,
            "target_chr": "os9",
            "target_start": 274373,
            "target_end": 289908
        },
        {
            "source_chr": "os5",
            "source_start": 337479,
            "source_end": 413516,
            "target_chr": "os8",
            "target_start": 456526,
            "target_end": 459280
        },
        {
            "source_chr": "os9",
            "source_start": 68195,
            "source_end": 140991,
            "target_chr": "os4",
            "target_start": 394943,
            "target_end": 412270
        },
        {
            "source_chr": "os14",
            "source_start": 184401,
            "source_end": 411744,
            "target_chr": "os9",
            "target_start": 483605,
            "target_end": 588244
        },
        {
            "source_chr": "os8",
            "source_start": 679191,
            "source_end": 686885,
            "target_chr": "os13",
            "target_start": 744636,
            "target_end": 1016986
        },
        {
            "source_chr": "os12",
            "source_start": 539417,
            "source_end": 728038,
            "target_chr": "os5",
            "target_start": 659441,
            "target_end": 782888
        },
        {
            "source_chr": "os5",
            "source_start": 217047,
            "source_end": 227875,
            "target_chr": "os8",
            "target_start": 187312,
            "target_end": 297956
        }
    ]
    const finalChords = (chordData.map((chord) => ({
            source_chromosome: chord.source_chr,
            source_start: chord.source_start,
            source_end: chord.source_end,
            target_chromosome: chord.target_chr,
            target_start: chord.target_start,
            target_end: chord.target_end,
    })));
    

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
    };

    const handleChordConfigUpdate = (updatedConfig: Partial<ChordConfig>) => {
        setChordConfig((prevConfig) => ({ ...prevConfig, ...updatedConfig }));
    }

    const handleGlobalConfigUpdate = (updatedConfig: Partial<ChordConfig>) => {
        setGlobalConfig((prevConfig) => ({ ...prevConfig, ...updatedConfig }));
    }

    return (
        <div>
            <FileUpload onFileUpload={handleFileUpload} />
            <div style={{ display: "flex", flexDirection: "row", height: "100vh" }}>
                <div style={{ flex: 5, display: "flex", alignItems: "center", justifyContent: "center" }}>
                    <div ref={canvasRef} style={{ border: "1px solid black" }}>
                        <Segment data={{ segments, config: assemblyConfig, globalConfig: globalConfig, divRef: canvasRef }} />
                        <Chords data={{ chords: finalChords, config: chordConfig, globalConfig: globalConfig, divRef: canvasRef }} />
                    </div>
                </div>

                <div style={{ flex: 5, display: "flex", alignItems: "center", justifyContent: "center" }}>
                    {segments.length > 0 && (
                        <Form handleAssemblyConfigUpdate={handleAssemblyConfigUpdate} defaultAssemblyConfig={assemblyConfig} 
                            handleChordConfigUpdate={handleChordConfigUpdate} defaultChordConfig={chordConfig}
                            handleGlobalConfigUpdate={handleGlobalConfigUpdate} defaultGlobalConfig={globalConfig}
                        />
                    )}
                </div>
            </div>
        </div>
    );
};

export default Circos;
