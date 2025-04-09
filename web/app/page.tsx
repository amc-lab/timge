"use client";
import React, { useEffect } from "react";
import "@fontsource/roboto";
import Header from "@/components/Header";
import { v4 as uuidv4 } from "uuid";
import CircosView from "./circos/CircosView";
import { useState } from "react";
import TrackUploadForm from "@/components/TrackUploadForm";
import { State } from "./types/state";
import { ViewType } from "./types/viewTypes";
import MapView from "./map/MapView";
import LinearView from "./linear/LinearView";
import { Box } from "@mui/joy";
import { Hash } from "crypto";

interface FileEntry {
  name: string;
  data: any;
  trackType: string;
}

function getDefaultState(): State {
  const defaultUUID = uuidv4();
  return {
    title: 'Untitled',
    // dateCreated: new Date().toISOString(),
    // dateModified: new Date().toISOString(),
    views: [],
    isUserLoggedIn: false,
    UUID: defaultUUID,
    // dataFilesRootDir: `data/${defaultUUID}/`,
    dataFiles: [],
  };
}

export default function Page() {

  const fileFormatMapping = {
    "fasta": "karyotype",
    "bed": "line",
    "bedgraph": "line",
    "bedpe": "link",
    "fa": "karyotype",
  }

  const STATE_KEY = "timge-state";

  function loadState() {
    if (typeof window === "undefined") return getDefaultState();
    const savedState = localStorage.getItem(STATE_KEY);
    if (savedState) {
      try {
        return JSON.parse(savedState);
      }
      catch (error) {
        console.error("Error parsing state from localStorage:", error);
        return getDefaultState();
      }
    }
    return getDefaultState();
  }

  function saveState(newState: any) {
    if (typeof window !== "undefined") {
      localStorage.setItem(STATE_KEY, JSON.stringify(newState));
    }
  }

  function exportState() {
    const state = loadState();
    const blob = new Blob([JSON.stringify(state)], { type: "application/json" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "timge-state.json";
    a.click();
    URL.revokeObjectURL(url);
  }

  function deleteTrack(trackName: string) {
    setTracks(tracks.filter((track) => track.name !== trackName));
    setSpaceState({
      ...spaceState,
      dataFiles: spaceState.dataFiles.filter((file) => file !== trackName),
    });
    const host = process.env.NEXT_PUBLIC_DJANGO_HOST;
    fetch(`${host}/api/timge/delete_track/?uuid=${spaceState.UUID}&track_name=${trackName}`, {
      method: "DELETE",
    })
    .then((response) => response.json())
    .then((data) => {
      if (data.status === "success") {
        console.log("Track deleted successfully");
      } else {
        console.error("Error deleting track:", data.message);
      }
    })
  }

  const [files, setFiles] = useState<FileEntry[]>([]);
  const [tracks, setTracks] = useState([]);
  const [genomeFormOpen, setGenomeFormOpen] = useState(false);
  const [spaceState, setSpaceState] = useState<State>(() => loadState());
  // Mapping from view ID to linked view ID
  const [connections, setConnections] = useState<Map<string, string[]>>(new Map());
  // Maintains dependency properties for each view
  const [dependencies, setDependencies] = useState<Map<string, any>>(new Map());
  const [createdViews, setCreatedViews] = useState<Set<string>>(new Set());

  // Update the space state in local storage whenever it changes
  useEffect(() => {
    saveState(spaceState);
  }, [spaceState]);

  // Get track files from the backend on initial load
  useEffect(() => {
    const getTrackFiles = () => {
      const host = process.env.NEXT_PUBLIC_DJANGO_HOST;
      fetch(`${host}/api/timge/get_tracks/?uuid=${spaceState.UUID}`, {
        method: "GET",
      })
        .then((response) => response.json())
        .then((data) => {
          if (data.status === "success") {
            const trackFiles = data.track_files.map((file) => ({
              name: file.name,
              data: file.content,
              type: file.type,
              size: file.size,
              trackType: fileFormatMapping[file.name.split('.').pop()],
            }));
  
            setFiles(trackFiles);
  
          } else {
            console.error("Error fetching tracks:", data.message);
          }
        })
        .catch((err) => {
          console.error("Request failed:", err);
        });
    };
  
    getTrackFiles();
  }, []);

  // Update the tracks whenever the files change
  useEffect(() => {
    let circosFiles = [];

    const formData = new FormData();
    formData.append("track_types", JSON.stringify(files.map((file) => fileFormatMapping[file.name.split(".").pop()])));
    files.forEach((file) => {
      formData.append("data_files", new Blob([file.data]), file.name);
    });

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
          name: files[i].name,
          trackType: fileFormatMapping[files[i].name.split(".").pop()],
        });
      }      
    })
    .then(() => {
      setTracks([...tracks, ...circosFiles]);
    });
  }
  , [files]);

  useEffect(() => {
    console.log("Connections:", connections);
    console.log("Dependencies:", dependencies);
    console.log("Created Views:", createdViews);
  }, [connections, dependencies, createdViews]);

  const addLinearGenomeView = () => {
    const id = `view-${spaceState.views.length + 1}`;
    const newView = {
      viewType: "LinearGenomeView",
    }
    setSpaceState({
      ...spaceState,
      views: [
        ...spaceState.views,
        {
          id: `view-${spaceState.views.length + 1}`,
          type: ViewType.Linear,
          title: `Linear Genome View ${spaceState.views.length + 1}`,
          description: "Standard linear genome view",
          uuid: spaceState.UUID,
          config: {
            isMinimised: false,
          }
        },
      ],
    });
    setCreatedViews((prev) => new Set(prev).add(id));
    setDependencies((prev) => new Map(prev).set(id, {}));
    setConnections((prev) => new Map(prev).set(id, []));
  }

  const addCircosView = () => {
    const id = `view-${spaceState.views.length + 1}`;
    setSpaceState({
      ...spaceState,
      views: [
        ...spaceState.views,
        {
          id: id,
          type: ViewType.Circos,
          title: `Circos View ${spaceState.views.length + 1}`,
          description: "Circular genome view",
          visible_tracks: [],
          uuid: spaceState.UUID,
          config: {
            isMinimised: false,
          },
        },
      ],
    });
    setDependencies((prev) => new Map(prev).set(id, {}));
    setConnections((prev) => new Map(prev).set(id, []));
    setCreatedViews((prev) => new Set(prev).add(id));
  }

  const addMapView = () => {
    const id = `view-${spaceState.views.length + 1}`;
    setSpaceState({
      ...spaceState,
      views: [
        ...spaceState.views,
        {
          id: id,
          type: ViewType.Map,
          title: `Map View ${spaceState.views.length + 1}`,
          description: "Map view",
          uuid: spaceState.UUID,
          config: {
            reference: "",
            track: "",
            segmentA: "",
            segmentB: "",
            resolution: 5,
            isMinimised: false,
          },
        },
      ],
    });
    setDependencies((prev) => new Map(prev).set(id, {}));
    setConnections((prev) => new Map(prev).set(id, []));
    setCreatedViews((prev) => new Set(prev).add(id));
  }

  const addConnection = (viewId: string, linkedViewId: string) => {
    setConnections((prev) => {
      const newConnections = new Map(prev);
      newConnections.set(viewId, [...(newConnections.get(viewId) || []), linkedViewId]);
      return newConnections;
    });
  }

  const removeConnection = (viewId: string) => {
    const dependents = connections.get(viewId);
    if (dependents) {
      dependents.forEach((dependentId) => {
        setDependencies((prev) => {
          const newDependencies = new Map(prev);
          newDependencies.set(dependentId, {});
          return newDependencies;
        }
        );
      }
      );
    }

    setConnections((prev) => {
      const newConnections = new Map(prev);
      newConnections.delete(viewId);
      return newConnections;
    });
  }

  const addCustomMapView = (mapConfig) => {
    setSpaceState({
      ...spaceState,
      views: [
        ...spaceState.views,
        {
          id: `view-${spaceState.views.length + 1}`,
          type: ViewType.Map,
          title: `Map View ${spaceState.views.length + 1}`,
          description: "Map view",
          uuid: spaceState.UUID,
          config: mapConfig,
        },
      ],
    });
  }

  const crossViewActionHandler = (action: string, data: any) => {
    if (action === "generate_heatmap") {
      if (!data || !data.track || !data.reference || !data.segmentA || !data.segmentB || !data.resolution) {
        console.error("Invalid or missing data for generate_heatmap action");
        return;
      }
      console.log("Generating heatmap with data:", data);
      addCustomMapView({
        reference: data.reference,
        track: data.track,
        segmentA: data.segmentA,
        segmentB: data.segmentB,
        resolution: data.resolution,
      });
    }
    else if (action === "delete_view") {
      const viewId = data.viewId;
      if (viewId) {
        setSpaceState((prevState) => ({
          ...prevState,
          views: prevState.views.filter((view) => view.id !== viewId),
        }));
        setConnections((prev) => {
          const newConnections = new Map(prev);
          newConnections.delete(viewId);
          return newConnections;
        });
        setDependencies((prev) => {
          const newDependencies = new Map(prev);
          newDependencies.delete(viewId);
          return newDependencies;
        });
        setCreatedViews((prev) => {
          const newCreatedViews = new Set(prev);
          newCreatedViews.delete(viewId);
          return newCreatedViews;
        }
        );
      }
    }
    else if (action === "add_connection") {
      const { source, target } = data;
      if (source && target) {
        addConnection(source, target);
      }
    }
    else if (action === "remove_connection") {
      const { viewId } = data;
      if (viewId) {
        removeConnection(viewId);
      }
    }
    else if (action === "propagate_dependencies") {
      const { viewId, dependencies } = data;
      console.log("Propagating dependencies for view:", viewId, "with dependencies:", dependencies);
      const dependents = connections.get(viewId);
      console.log("Dependents:", dependents);
      if (dependents) {
        dependents.forEach((dependentId) => {
          console.log("Setting dependencies for dependent view:", dependentId);
          setDependencies((prev) => {
            const newDependencies = new Map(prev);
            newDependencies.set(dependentId, dependencies);
            return newDependencies;
          });
        });
      }
    }
    else {
      console.log("Unknown action:", action);
    }
  }

  // When importing a state, set the space state to the imported state and fetch the tracks from the backend
  const importState = () => {
    const fileInput = document.createElement("input");
    fileInput.type = "file";
    fileInput.accept = ".json";
    fileInput.onchange = (event) => {
      const file = (event.target as HTMLInputElement).files?.[0];
      if (file) {
        const reader = new FileReader();
        reader.onload = (e) => {
          const contents = e.target?.result;
          if (contents) {
            const state = JSON.parse(contents as string);
            console.log("Loaded state:", state);
            setSpaceState(state);
            localStorage.setItem(STATE_KEY, JSON.stringify(state));
            // fetch the tracks from the backend
            const host = process.env.NEXT_PUBLIC_DJANGO_HOST;
            fetch(`${host}/api/timge/get_tracks/?uuid=${state.UUID}`, {
              method: "GET",
            })
            .then((response) => response.json())
            .then((data) => {
              if (data.status === "success") {
                const trackFiles = data.track_files.map((file) => ({
                  name: file.name,
                  data: file.content,
                  type: file.type,
                  size: file.size,
                  trackType: fileFormatMapping[file.name.split('.').pop()],
                }));
                setFiles(trackFiles);
              } else {
                console.error("Error fetching tracks:", data.message);
              }
            }
            )
          }
        };
        reader.readAsText(file);
      }
    };
    fileInput.click();
    const state = loadState();
    setSpaceState(state);
  }

  const importTracks = () => {
    setGenomeFormOpen(true);
  }

  const updateViewState = (index: number, updatedConfig: any) => {
    setSpaceState((prevState) => {
      const updatedViews = [...prevState.views];
      updatedViews[index] = {
        ...updatedViews[index],
        ...updatedConfig,
      };
      return {
        ...prevState,
        views: updatedViews,
      };
    });
  }

  return <>
    <Header 
      addLinearGenomeView={addLinearGenomeView}
      addCircosView={addCircosView}
      addMapView={addMapView}
      importTracks={importTracks}
      importState={importState}
      exportState={exportState}
      resetState={() => {
        // make a request to the backend to delete the files
        const host = process.env.NEXT_PUBLIC_DJANGO_HOST;
        fetch(`${host}/api/timge/delete_tracks/?uuid=${spaceState.UUID}`, {
          method: "DELETE",
        })
        .then((response) => response.json())
        .then((data) => {
          if (data.status === "success") {
            console.log("Tracks deleted successfully");
          } else {
            console.error("Error deleting tracks:", data.message);
          }
        })
        localStorage.removeItem(STATE_KEY);
        setTracks([]);
        setSpaceState(getDefaultState());
      }}
    />
    
    <TrackUploadForm 
      isOpen={genomeFormOpen} 
      onClose={() => setGenomeFormOpen(false)} 
      tracks={tracks}
      onDeleteTrack={deleteTrack}
      onTrackUpload={(data_files: File[]) => {
        const newTracks = data_files.map((file) => ({
          name: file.name,
          data: file,
          trackType: fileFormatMapping[file.name.split(".").pop()],
        }));
        setFiles(newTracks);
        setSpaceState({
          ...spaceState,
          dataFiles: [
            ...spaceState.dataFiles,
            ...newTracks.map((track) => track.name),
          ],
        });
        const formData = new FormData();
        newTracks.forEach((track) => {
          formData.append("track_files", track.data);
        }
        );
        formData.append("uuid", spaceState.UUID);
        const host = process.env.NEXT_PUBLIC_DJANGO_HOST;
        fetch(`${host}/api/timge/upload_tracks/`, {
          method: "POST",
          body: formData,
        })
        .then((response) => response.json())
        .then((data) => {
          console.log("Files uploaded successfully", data);
        }
        )
        .catch((error) => {
          console.error("Error uploading files", error);
        }
        );
      }}
    />
    <Box
      sx={{
        display: "flex",
        flexDirection: "row",
        flexWrap: "wrap",
        alignItems: "flex-start",
      }}
      >
      {
        spaceState.views.map((view, index) => {
          if (view.type === "linear") {
            return <LinearView 
                      key={index} 
                      trackFiles={tracks} 
                      viewConfig={view} 
                      handleViewUpdate={updateViewState} 
                      index={index} 
                      crossViewActionHandler={crossViewActionHandler} 
                      dependencies={dependencies.get(view.id)}
                      addConnection={addConnection}
                      removeConnection={removeConnection}
                      createdViews={createdViews}
                    />  
          }
          else if (view.type === "circos") {
            return <CircosView 
                      key={index}
                      trackFiles={tracks}
                      viewConfig={view}
                      handleViewUpdate={updateViewState} 
                      index={index} 
                      crossViewActionHandler={crossViewActionHandler} 
                      dependencies={dependencies.get(view.id)}
                      addConnection={addConnection}
                      removeConnection={removeConnection}
                      createdViews={createdViews}
                    />
          }
          else if (view.type === "map") {
            return <MapView 
                      key={index} 
                      trackFiles={tracks} 
                      viewConfig={view}
                      handleViewUpdate={updateViewState} 
                      index={index} 
                      crossViewActionHandler={crossViewActionHandler} 
                      dependencies={dependencies.get(view.id)}
                      addConnection={addConnection}
                      removeConnection={removeConnection}
                      createdViews={createdViews}
                    />
          }
        })
      }
    </Box>
    </>;

}