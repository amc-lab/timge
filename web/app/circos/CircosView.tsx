"use client"
import { Box, Button, Card, IconButton, Input, Option, Select, Slider, TextField, Typography } from "@mui/joy";
import { Track, TrackType } from "./config/track";
import { useState } from "react";
import Tracks from "./tracks";
import { defaultAnnotationConfig, defaultAssemblyConfig, defaultChordConfig, defaultGlobalConfig, defaultLineConfig, defaultHighlightConfig } from "./config/defaultConfigs";
import { useRef, useEffect } from "react";
import MenuIcon from "@mui/icons-material/Menu";
import * as d3 from "d3";
import ParentView from "@/components/ParentView";
import TrackSelector from "./components/TrackSelector";
import { View } from "@/store/features/views/types";
import { useAppDispatch, useAppSelector } from "@/store/hooks";
import { setConnection, setDependency } from "@/store/features/space/spaceSlice";

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
  const [connectedViews, setConnectedViews] = useState<string[]>(() => {
    const direct = space.connections[props.viewConfig.uuid] || [];
    const reverse = Object.entries(space.connections)
      .filter(([key, value]) =>
        Array.isArray(value) && value.includes(props.viewConfig.uuid)
      )
      .map(([key]) => key);
    return Array.from(new Set([...direct, ...reverse]));
  });
  const [minFilterScore, setMinFilterScore] = useState(0);

  const fileFormatMapping = {
    "fasta": "karyotype",
    "bed": "annotation",
    "bedgraph": "line",
    "bedpe": "link",
    "fa": "karyotype",
  }
  
  const [tracks, setTracks] = useState<Track[]>([]);
  const [selectedTracks, setSelectedTracks] = useState<Track[]>([]);

  const shouldDisplayFilterScore = (tracks: Track[]) => {
    return tracks.some((track) => track.trackType === TrackType.Chord);
  };
  const getMaxFilterScore = (tracks: Track[]) => {
    const chordTracks = tracks.filter((track) => track.trackType === TrackType.Chord);
    if (chordTracks.length > 0) {
      return chordTracks[0].data.chords.reduce((max, chord) => Math.max(max, chord.score), 0) + 1;
    }
    return 1000;
  };

  const generate_circos_files = (files: any) => {
    console.log("Generating circos files", files);
    const FILE_HOST = process.env.NEXT_PUBLIC_FILE_HOST;
    let circosFiles = [];

    const formData = new FormData();
    formData.append("track_types", JSON.stringify(files.map((file) => fileFormatMapping[file.split(".").pop()])));
    formData.append("track_paths", JSON.stringify(files.map((file) => FILE_HOST + space.uuid + "/" + file)));
    
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
    let updatedTracks: Track[] = [
      {
        trackType: TrackType.Highlight,
        config: defaultHighlightConfig,
        data: {
          globalConfig,
          divRef: canvasRef,
        },
        name: "Highlight",
      }
    ];
    files.forEach((trackFile) => {
        if (trackFile.name.endsWith(".bedgraph")) {
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
        else if (trackFile.name.endsWith(".bed")) {
            console.log("Annotation track", trackFile, globalConfig, canvasRef);
            updatedTracks.push({
                trackType: TrackType.Annotation,
                config: defaultAnnotationConfig,
                data: {
                    annotations: trackFile.data,
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
            "Export SVG": () => {
              const svg = d3.select(canvasRef.current).select("svg");
              const serializer = new XMLSerializer();
              const svgString = serializer.serializeToString(svg.node());
              const blob = new Blob([svgString], { type: "image/svg+xml" });
              const url = URL.createObjectURL(blob);
              const a = document.createElement("a");
              a.href = url;
              a.download = props.viewConfig.title + ".svg";
              document.body.appendChild(a);
              a.click();
              document.body.removeChild(a);
              URL.revokeObjectURL(url);
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
              onClose={handleTrackSelectorClose}
              onConfirm={(selectedTracks) => {
                props.handleViewUpdate(props.index, {
                  ...props.viewConfig,
                  visible_tracks: selectedTracks.map((track) => track.path),
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
                <Slider
                  value={minFilterScore}
                  min={0}
                  max={getMaxFilterScore(selectedTracks)}
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
                    Link Transparency:
                  </Typography>
                  <Slider
                    min={0}
                    max={1}
                    step={0.1}
                    value={globalConfig.linkUnselectedOpacity}
                    valueLabelDisplay="auto"
                    sx={{
                      width: "10%",
                      marginLeft: "10px",
                    }}
                    onChange={(event, newValue) => {
                      setGlobalConfig({
                        ...globalConfig,
                        linkUnselectedOpacity: newValue as number,
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
                      // Find added and removed connections
                      const prev = connectedViews;
                      const added = value.filter((v) => !prev.includes(v));
                      const removed = prev.filter((v) => !value.includes(v));

                      setConnectedViews(value);

                      // Add new connections
                      for (let i = 0; i < added.length; i++) {
                        const source = props.viewConfig.uuid;
                        const target = added[i];
                        const targetView = space.views.find((view) => view.uuid === target);
                        if (targetView?.type === "linear") {
                          dispatch(setConnection({
                            key: target,
                            value: [...(space.connections[target] || []), source]
                          }));
                        } else {
                          dispatch(setConnection({
                            key: source,
                            value: [...(space.connections[source] || []), target]
                          }));
                        }
                      }

                      // Remove deselected connections
                      for (let i = 0; i < removed.length; i++) {
                        const source = props.viewConfig.uuid;
                        const target = removed[i];
                        const targetView = space.views.find((view) => view.uuid === target);
                        if (targetView?.type === "linear") {
                          // Remove source from target's connections
                          const updated = (space.connections[target] || []).filter((uuid) => uuid !== source);
                          dispatch(setConnection({
                            key: target,
                            value: updated
                          }));
                          dispatch(setDependency({
                            key: source,
                            value: null
                          }));
                        } else {
                          // Remove target from source's connections
                          const updated = (space.connections[source] || []).filter((uuid) => uuid !== target);
                          dispatch(setConnection({
                            key: source,
                            value: updated
                          }));
                        }
                      }
                    }}
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
                <Tracks 
                  id={props.viewConfig.uuid} 
                  tracks={selectedTracks} 
                  crossViewActionHandler={props.crossViewActionHandler} 
                  globalConfig={globalConfig} 
                  dependencies={props.dependencies}
                  />
              </div>
              </Box>
            </Box>
          )}
        </ParentView>
      );
}

export default CircosView;