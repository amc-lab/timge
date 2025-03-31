"use client"
import { Box, IconButton } from "@mui/joy";
import { Track, TrackType } from "./config/track";
import { useState } from "react";
import Tracks from "./tracks";
import { defaultAssemblyConfig, defaultChordConfig, defaultGlobalConfig, defaultLineConfig } from "./config/defaultConfigs";
import { useRef, useEffect } from "react";
import MenuIcon from "@mui/icons-material/Menu";
import * as d3 from "d3";

interface CircosViewProps {
    trackFiles: any[];
}

const CircosView = (props: CircosViewProps) => {
   const canvasRef = useRef<HTMLDivElement>(null);
   const [globalConfig, setGlobalConfig] = useState(defaultGlobalConfig);

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
                    }
                });
            } else if (trackFile.name.endsWith(".fa")) {
                updatedTracks.push({
                    trackType: TrackType.Karyotype,
                    config: defaultAssemblyConfig,
                    data: {
                        segments: trackFile.data,
                        globalConfig,
                        divRef: canvasRef,
                    }
                });
            } else if (trackFile.name.endsWith(".bedpe")) {
                updatedTracks.push({
                    trackType: TrackType.Chord,
                    config: defaultChordConfig,
                    data: {
                        chords: trackFile.data,
                        globalConfig,
                        divRef: canvasRef,
                    }
                });
            }
        });
        setTracks(updatedTracks);
    }
    , [props.trackFiles]);
    
    return (
        <Box
            sx={{
                display: "flex",
                justifyContent: "center",
                alignItems: "center",
                width: "calc(100% - 5px)",
                borderRadius: "3px",
                margin: "2.5px",
                flexDirection: "column",
                backgroundColor: "white",
                border: "4px solid darkblue",
            }}
        >
            <Box
                sx={{
                    display: "flex",
                    // justifyContent: "center",
                    alignItems: "center",
                    width: "100%",
                    height: "2em",
                    backgroundColor: "darkblue",
                }}
            >
            <IconButton sx={{ color: "white", height: "2em", "&:hover": { background: "none", color: "white" } }}>
                <MenuIcon />
            </IconButton>
            </Box>
            <Box
                sx={{
                    backgroundColor: "white",
                    borderRadius: "3px",
                    width: "650px",
                    height: "650px",  
                    padding: "5px"  
                }}
            >
            <div ref={canvasRef} style={{ width: "100%", height: "100%" }}>
                <Tracks tracks={tracks} />
              </div>
            </Box>
        </Box>
    )
}

export default CircosView;
