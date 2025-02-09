"use client";
import {
  type Assembly,
  type Chord,
  type AssemblyConfig,
  type ChordConfig,
  type GlobalConfig,
} from "../types/genomes";
import { useRef, useEffect, useState } from "react";
import * as d3 from "d3";
import Form from "../../components/Form";
import Tracks from "./tracks";
import { Track, TrackType } from "./config/track";
import {
  defaultAssemblyConfig,
  defaultBarConfig,
  defaultChordConfig,
  defaultGlobalConfig,
} from "./config/defaultConfigs";
import FileUploadForm from "../../components/FileUploadForm";
import { Box, Button } from "@mui/joy";

interface CircosProps {
  params: Promise<{
    segments: Array<Assembly>;
    chords: Array<Chord>;
  }>;
}

export default function Circos({ params }: CircosProps) {
  const canvasRef = useRef<HTMLDivElement>(null);
  const [segments, setSegments] = useState<Assembly[]>([]);
  const [globalConfig, setGlobalConfig] = useState(defaultGlobalConfig);
  const [tracks, setTracks] = useState<Track[]>([]);

  useEffect(() => {
    params.then(({ segments, chords }) => {
      setSegments(segments);
      // Handle chords if needed
    });

    if (canvasRef.current) {
      const svg = d3.select(canvasRef.current).select("svg");

      if (svg.empty()) {
        d3.select(canvasRef.current)
          .append("svg")
          .attr("width", globalConfig.canvasWidth)
          .attr("height", globalConfig.canvasHeight);
      }
    }
  }, [params, globalConfig]);

  const handleGlobalConfigUpdate = (updatedConfig: Partial<GlobalConfig>) => {
    setGlobalConfig((prevConfig) => ({ ...prevConfig, ...updatedConfig }));
  };

  const handleTrackConfigUpdate = (index: number, updatedConfig: any) => {
    setTracks((prevTracks) => {
      const newTracks = [...prevTracks];
      newTracks[index] = { ...newTracks[index], config: updatedConfig };
      return newTracks;
    });
  };

  const onSubmit = (files) => {
    const newTracks = files.map((file) => {
      switch (file.trackType) {
        case "karyotype":
          setSegments(file.data);
          return {
            trackType: TrackType.Karyotype,
            data: { segments: file.data, globalConfig, divRef: canvasRef },
            config: defaultAssemblyConfig,
          };
        case "bar":
          return {
            trackType: TrackType.Bar,
            data: { bars: file.data, globalConfig, divRef: canvasRef },
            config: defaultBarConfig,
          };
        case "link":
          return {
            trackType: TrackType.Chord,
            data: { chords: file.data, globalConfig, divRef: canvasRef },
            config: defaultChordConfig,
          };
        default:
          throw new Error("Invalid track type");
      }
    });

    setTracks(newTracks);
  };

  const downloadSVGAsPNG = () => {
    const d3ToPng = require('d3-svg-to-png');
    d3ToPng(canvasRef.current, 'output', {
      scale: 1,
      format: 'png',
      quality: 1,
      download: true,
      ignore: '.ignored',
      background: 'white'
    });
}

  return (
    <div>
      <FileUploadForm onSubmit={onSubmit} />
      <div style={{ display: "flex", flexDirection: "row", height: "100vh" }}>
        <div
          style={{
            flex: 5,
            display: "flex",
            alignItems: "center",
            justifyContent: "center",
          }}
        >
          {
            segments && segments.length > 0 && (
              <div ref={canvasRef} style={{ border: "1px solid black" }}>
                <Tracks tracks={tracks} />
              </div>
            )
          }
          
        </div>

        <div
          style={{
            flex: 5,
            display: "flex",
            alignItems: "center",
            justifyContent: "center",
          }}
        >
          {segments && segments.length > 0 && (
          <Box>
            <Form
                tracks={tracks}
                handleTrackConfigUpdate={handleTrackConfigUpdate}
                handleGlobalConfigUpdate={handleGlobalConfigUpdate}
              />
          <Button onClick={downloadSVGAsPNG} style={{ margin: "10px", padding: "12px", fontSize: "12px" }}>
            Download PNG
          </Button>
          </Box>
          )}
        </div>
      </div>
    </div>
  );
}
