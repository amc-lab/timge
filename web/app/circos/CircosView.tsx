"use client"
import { Box, Button, Card, IconButton, Input, Option, Select, Slider, TextField, Typography } from "@mui/joy";
import { Track, TrackType } from "./config/track";
import { useState } from "react";
import Tracks from "./tracks";
import { defaultAssemblyConfig, defaultChordConfig, defaultGlobalConfig, defaultLineConfig } from "./config/defaultConfigs";
import { useRef, useEffect } from "react";
import MenuIcon from "@mui/icons-material/Menu";
import * as d3 from "d3";
import ParentView from "@/components/ParentView";
import TrackSelector from "./components/TrackSelector";
import { View } from "@/store/features/views/types";
import { useAppDispatch, useAppSelector } from "@/store/hooks";
import { setConnection } from "@/store/features/space/spaceSlice";

interface CircosViewProps {
    viewConfig: View;
    handleViewUpdate: (index, viewState: View) => void;
    crossViewActionHandler?: any;
    index: number;
    dependencies?: any;
    addConnection?: any;
    removeConnection?: any;
    connections?: string[];
    files: any[];
}

const CircosView = (props: CircosViewProps) => {
  const space = useAppSelector((state) => state.space);
  const dispatch = useAppDispatch();

  const canvasRef = useRef<HTMLDivElement>(null);
  const [globalConfig, setGlobalConfig] = useState(defaultGlobalConfig);
  const [connectedViews, setConnectedViews] = useState<string[]>([]);
  const [minFilterScore, setMinFilterScore] = useState(0);

  const fileFormatMapping = {
    "fasta": "karyotype",
    "bed": "line",
    "bedgraph": "line",
    "bedpe": "link",
    "fa": "karyotype",
  }
  
  const [tracks, setTracks] = useState<Track[]>([]);
  const [selectedTracks, setSelectedTracks] = useState<Track[]>([]);

  const generate_circos_files = (files: any) => {
    console.log("Generating circos files", files);
    let circosFiles = [];

    const formData = new FormData();
    formData.append("track_types", JSON.stringify(files.map((file) => fileFormatMapping[file.split(".").pop()])));
    formData.append("track_paths", JSON.stringify(files.map((file) => "/Users/mithun/Documents/Imperial/Year 4/FYP/uploaded_data/" + space.uuid + "/" + file)));

    const host = process.env.NEXT_PUBLIC_DJANGO_HOST;
    fetch(`${host}/api/multilift/circos/`, {
      method: "POST",
      body: formData,
    })
    .then((response) => response.json())
    .then((data) => {
      for (let i = 0; i < data.length; i++) {
        circosFiles.push({
          data: data[i],
          name: files[i],
          trackType: fileFormatMapping[files[i].split(".").pop()],
        });
      }      
    })
    .then(() => {
      const formattedCircosFiles = helper(circosFiles);
      setSelectedTracks([...tracks, ...formattedCircosFiles]);
    });
  }

  useEffect(() => {
    console.log("CircosView props", props.viewConfig);
    if (props.viewConfig.visible_tracks.length > 0) {
      generate_circos_files(props.viewConfig.visible_tracks);
    }
  }, [props.viewConfig.visible_tracks]);

  // useEffect(() => {
  //   if (props.files.length > 0) {
  //     generate_circos_files(props.files);
  //   }
  // }
  // , [props.files]);

  //  const maxScore = d3.max(props.trackFiles, (d) => {
  //      if (d.data) {
  //          return d3.max(d.data, (d) => d.score);
  //      }
  //      return 0;
  //  });

  //   console.log("CircosView props", props.trackFiles);

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

  const helper = (files) => {
    let updatedTracks: Track[] = [];
    files.forEach((trackFile) => {
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
        else {
          return;
        }
    });
    return updatedTracks;
  }

    const [isTrackSelectorOpen, setIsTrackSelectorOpen] = useState(false);
    const handleTrackSelectorClose = () => {
        setIsTrackSelectorOpen(false);
    };

    // useEffect(() => {
    //     // const selectedTracks = tracks.filter((track) => {
    //     //     return props.viewConfig.visible_tracks.includes(track.name);
    //     // });
    //     let selectedTracks: Track[] = [];
    //     for (let i = 0; i < props.viewConfig.visible_tracks.length; i++) {
    //         if (tracks.filter((track) => track.name === props.viewConfig.visible_tracks[i]).length > 0) {
    //             selectedTracks.push(tracks.filter((track) => track.name === props.viewConfig.visible_tracks[i])[0]);
    //         }
    //     }
    //     setSelectedTracks(selectedTracks);
    // }, [tracks]);

    return (
        <ParentView
          viewConfig={props.viewConfig}
          index={props.index}
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
              files={props.files}
              onClose={handleTrackSelectorClose}
              onConfirm={(selectedTracks) => {
                props.handleViewUpdate(props.index, {
                  ...props.viewConfig,
                  visible_tracks: selectedTracks.map((track) => track.name),
                });
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
                {/* <Input
                  type="number"
                  value={minFilterScore}
                  onChange={(event) => {
                    setMinFilterScore(parseInt(event.target.value));
                    setSelectedTracks((prevTracks) => {
                      return prevTracks.map((track) => {
                        if (track.trackType === TrackType.Chord) {
                          return {
                            ...track,
                            config: {
                              ...track.config,
                              minFilterScore: parseInt(event.target.value),
                            },
                          };
                        }
                        return track;
                      });
                    });
                  }}
                  sx={{
                    width: "5%",
                    marginLeft: "10px",
                  }}
                  ></Input> */}
                <Slider
                  // defaultValue={[20, maxScore]}
                  value={minFilterScore}
                  min={0}
                  max={1000}
                  step={1}
                  valueLabelDisplay="auto"
                  sx={{
                    width: "10%",
                  }}
                  onChange={(event, newValue) => {
                    setMinFilterScore(newValue as number);
                    setSelectedTracks((prevTracks) => {
                      return prevTracks.map((track) => {
                        if (track.trackType === TrackType.Chord) {
                          return {
                            ...track,
                            config: {
                              ...track.config,
                              minFilterScore: newValue as number,
                              // maxFilterScore: newValue[1],
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
                        // if (props.crossViewActionHandler) {
                          // console.log("Adding connection", value[i]);
                          // props.crossViewActionHandler("add_connection", {
                          //   source: props.viewConfig.uuid,
                          //   target: value[i],
                          // });
                          const source = props.viewConfig.uuid;
                          const target = value[i];
                          dispatch(setConnection(
                            { 
                              key: source, value: [...(space.connections[source] || []), target] 
                            }));
                        // }
                      }
                    }
                    }
      
                  >
                    {space.views &&
                      space.views.filter(
                        (view) => view.uuid !== props.viewConfig.uuid
                      )
                        .map((view) => {
                        return (
                          <Option
                            key={view.uuid}
                            value={view.uuid}
                          >
                            {view.title}
                          </Option>
                        );
                      }
                    )}
                  </Select>
              </Card>
              <Box
                sx={{
                  width: "650px",
                  paddingTop: "25px",
                  paddingBottom: "25px",
                }}
                >
              <div ref={canvasRef} style={{ width: "100%", height: "100%" }}>
                <Tracks id={props.viewConfig.uuid} tracks={selectedTracks} crossViewActionHandler={props.crossViewActionHandler} />
              </div>
              </Box>
            </Box>
          )}
        </ParentView>
      );
}

export default CircosView;