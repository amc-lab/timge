"use client"
import { Box, Button, Card, IconButton, Option, Select, Slider, Typography } from "@mui/joy";
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
    crossViewActionHandler?: any;
    index: number;
    dependencies?: any;
    addConnection?: any;
    removeConnection?: any;
    connections?: string[];
    createdViews: Set<any>;
}

const CircosView = (props: CircosViewProps) => {
   const canvasRef = useRef<HTMLDivElement>(null);
   const [globalConfig, setGlobalConfig] = useState(defaultGlobalConfig);
   const [connectedViews, setConnectedViews] = useState<string[]>([]);

   const maxScore = d3.max(props.trackFiles, (d) => {
       if (d.data) {
           return d3.max(d.data, (d) => d.score);
       }
       return 0;
   });

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
            } else if (trackFile.name.endsWith(".fa") || trackFile.name.endsWith(".fasta")) {
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

    return (
        <ParentView
          viewConfig={props.viewConfig}
          index={props.index}
          crossViewActionHandler={props.crossViewActionHandler}
          userActions={{
            "Select Tracks": () => {
              setIsTrackSelectorOpen(true);
            },
            "Download PNG": () => {
              const d3ToPng = require('d3-svg-to-png');
              d3ToPng(canvasRef.current, 'output', {
                scale: 1,
                format: 'png',
                quality: 1,
                download: true,
                ignore: '.ignored',
                background: 'white'
              });
            },
            "Clear": () => {
                props.handleViewUpdate(props.index, {
                    ...props.viewConfig,
                    visible_tracks: [],
                });
            },
            [props.viewConfig.config.isMinimised ? "Maximise" : "Minimise"]: () => {
              props.handleViewUpdate(props.index, {
                ...props.viewConfig,
                config: {
                  ...props.viewConfig.config,
                  isMinimised: !props.viewConfig.config.isMinimised,
                },
              });
            }
          }}
        >
          {isTrackSelectorOpen ? (
            <TrackSelector
              tracks={tracks}
              trackFiles={props.trackFiles}
              onClose={handleTrackSelectorClose}
              onConfirm={(selectedTracks) => {
                props.handleViewUpdate(props.index, {
                  ...props.viewConfig,
                  visible_tracks: selectedTracks.map((track) => track.name),
                });
                setSelectedTracks(selectedTracks);
                setIsTrackSelectorOpen(false);
              }}
            />
          ) : props.viewConfig.visible_tracks.length === 0 ? (
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
          ) : (
            <Box
              sx={{
                width: "100%",
                justifyContent: "center",
                display: "flex",
                flexDirection: "column",
                alignItems: "center",
              }}
            >
              <Card
                sx={{
                  width: "100%",
                  borderRadius: "none",
                  flexDirection: "row",
                  justifyContent: "center",
                  alignItems: "center",
                }}
              >
                <Typography
                  sx={{
                    fontSize: "1em",
                  }}
                >
                  Filter score:
                </Typography>
                <Slider
                  defaultValue={[0, maxScore]}
                  min={0}
                  max={maxScore}
                  step={1}
                  valueLabelDisplay="auto"
                  sx={{
                    width: "10%",
                  }}
                  onChange={(event, newValue) => {
                    setSelectedTracks((prevTracks) => {
                      return prevTracks.map((track) => {
                        if (track.trackType === TrackType.Chord) {
                          return {
                            ...track,
                            config: {
                              ...track.config,
                              minFilterScore: newValue[0],
                              maxFilterScore: newValue[1],
                            },
                          };
                        }
                        return track;
                      });
                    });
                  }}
                  ></Slider>
                  <Typography
                    sx={{
                      fontSize: "1em",
                      marginLeft: "10px",
                    }}
                  >
                    Connect to:
                  </Typography>
                  <Select
                    placeholder="Select views"
                    multiple
                    value={connectedViews}
                    onChange={(event, value) => {
                      setConnectedViews(value);
                      for (let i = 0; i < value.length; i++) {
                        if (props.crossViewActionHandler) {
                          props.crossViewActionHandler("add_connection", {
                            source: props.viewConfig.id,
                            target: value[i],
                          });
                        }
                      }
                    }
                    }
      
                  >
                    {props.createdViews &&
                      Array.from(props.createdViews).filter((view) => view.id !== props.viewConfig.id).map((view_id) => {
                        return (
                          <Option
                            key={view_id}
                            value={view_id}
                          >
                            {view_id}
                          </Option>
                        );
                      }
                    )}
                  </Select>
              </Card>
              <Box
                sx={{
                  width: "650px"
                }}
                >
              <div ref={canvasRef} style={{ width: "100%", height: "100%" }}>
                <Tracks id={props.viewConfig.id} tracks={selectedTracks} crossViewActionHandler={props.crossViewActionHandler} />
              </div>
              </Box>
            </Box>
          )}
        </ParentView>
      );
}

export default CircosView;