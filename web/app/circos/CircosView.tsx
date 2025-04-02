"use client"
import { Box, Button, IconButton } from "@mui/joy";
import { Track, TrackType } from "./config/track";
import { useState } from "react";
import Tracks from "./tracks";
import { defaultAssemblyConfig, defaultChordConfig, defaultGlobalConfig, defaultLineConfig } from "./config/defaultConfigs";
import { useRef, useEffect } from "react";
import MenuIcon from "@mui/icons-material/Menu";
import * as d3 from "d3";
import ParentView from "@/components/ParentView";
import TrackSelector from "./components/TrackSelector";
import { View } from "../types/state";

interface CircosViewProps {
    trackFiles: any[];
    viewConfig: View;
    handleViewUpdate: (index, viewState: View) => void;
    index: number;
}

const CircosView = (props: CircosViewProps) => {
   const canvasRef = useRef<HTMLDivElement>(null);
   const [globalConfig, setGlobalConfig] = useState(defaultGlobalConfig);

    console.log("CircosView props", props.trackFiles);

     useEffect(() => {
       if (canvasRef.current) {
         const svg = d3.select(canvasRef.current).select("svg");
   
         if (svg.empty()) {
           d3.select(canvasRef.current)
             .append("svg")
             .attr("width", "100%")
             .attr("height", "100%")
             .attr("viewBox", `0 0 ${globalConfig.canvasWidth} ${globalConfig.canvasHeight}`)
            .attr("preserveAspectRatio", "xMinYMin meet");
         }
        }
     }, [globalConfig]);

   const [tracks, setTracks] = useState<Track[]>([]);
   const [selectedTracks, setSelectedTracks] = useState<Track[]>([]);

    useEffect(() => {
        console.log("Setting track data", props.trackFiles);
        if (!props.trackFiles) {
            return;
        }
        let updatedTracks: Track[] = [];
        props.trackFiles.forEach((trackFile) => {
            if (trackFile.name.endsWith(".bed") || trackFile.name.endsWith(".bedgraph")) {
                updatedTracks.push({
                    trackType: TrackType.Line,
                    config: defaultLineConfig,
                    data: {
                        values: trackFile.data,
                        globalConfig,
                        divRef: canvasRef,
                    },
                    name: trackFile.name,
                });
            } else if (trackFile.name.endsWith(".fa")) {
                updatedTracks.push({
                    trackType: TrackType.Karyotype,
                    config: defaultAssemblyConfig,
                    data: {
                        segments: trackFile.data,
                        globalConfig,
                        divRef: canvasRef,
                    },
                    name: trackFile.name,
                });
            } else if (trackFile.name.endsWith(".bedpe")) {
                updatedTracks.push({
                    trackType: TrackType.Chord,
                    config: defaultChordConfig,
                    data: {
                        chords: trackFile.data,
                        globalConfig,
                        divRef: canvasRef,
                    },
                    name: trackFile.name,
                });
            }
        });
        setTracks(updatedTracks);
    }
    , [props.trackFiles]);

    const [isTrackSelectorOpen, setIsTrackSelectorOpen] = useState(false);
    const handleTrackSelectorClose = () => {
        setIsTrackSelectorOpen(false);
    };

    useEffect(() => {
        // const selectedTracks = tracks.filter((track) => {
        //     return props.viewConfig.visible_tracks.includes(track.name);
        // });
        let selectedTracks: Track[] = [];
        for (let i = 0; i < props.viewConfig.visible_tracks.length; i++) {
            if (tracks.filter((track) => track.name === props.viewConfig.visible_tracks[i]).length > 0) {
                selectedTracks.push(tracks.filter((track) => track.name === props.viewConfig.visible_tracks[i])[0]);
            }
        }
        setSelectedTracks(selectedTracks);
    }, [tracks]);

    if (props.viewConfig.visible_tracks.length === 0) {
        return (
            <ParentView
                viewConfig={props.viewConfig}
            >
                {
                    isTrackSelectorOpen && (
                        <TrackSelector
                            tracks={tracks}
                            trackFiles={props.trackFiles}
                            onClose={handleTrackSelectorClose}
                            onConfirm={(selectedTracks) => {
                                props.handleViewUpdate(props.index, {
                                    ...props.viewConfig,
                                    visible_tracks: selectedTracks.map((track) => track.name)
                                });
                                setSelectedTracks(selectedTracks);
                                setIsTrackSelectorOpen(false);
                            }}
                        />
                    )
                }
                <Box
                    sx={{
                        display: "flex",
                        justifyContent: "center",
                        alignItems: "center",
                        width: "100%",
                        height: "100%",
                        flexDirection: "column",
                        flexGrow: 1,
                        padding: 2,
                    }}
                >
                    <p>No Tracks Selected</p>
                    <Button
                        color="primary"
                        onClick={() => {
                            setIsTrackSelectorOpen(true);
                        }}
                        sx={{
                            marginTop: "10px",
                        }}
                    >
                        Select Tracks
                    </Button>
                </Box>
            </ParentView>
        );
    }
    
    return (
        <ParentView
            viewConfig={props.viewConfig}
        >
            <div ref={canvasRef} style={{ width: "100%", height: "100%" }}>
                <Tracks tracks={selectedTracks} />
              </div>
        </ParentView>
    )
}

export default CircosView;
