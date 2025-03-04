import { Box } from "@mui/joy";
import { Track, TrackType } from "./config/track";
import {type Assembly} from "../types/genomes";
import { useState } from "react";
import Tracks from "./tracks";

const CircosView = (props) => {
    
    return (
        <Box
            sx={{
                display: "flex",
                justifyContent: "center",
                alignItems: "center",
                width: "100%",
                borderRadius: "2em",
                margin: "5px",
            }}
        >
            <Box
                sx={{
                    display: "flex",
                    justifyContent: "center",
                    alignItems: "center",
                    width: "100%",
                    height: "2em",
                    backgroundColor: "lightblue",
                }}
            >
            </Box>
            <Tracks tracks={tracks} />
        </Box>
    )
}

export default CircosView;
