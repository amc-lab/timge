"use client";
export const dynamic = "force-dynamic";
import React, { useEffect, useState } from "react";
import "@fontsource/roboto";
import Header from "@/components/Header";
import CircosView from "./circos/CircosView";
import TrackUploadForm from "@/components/TrackUploadForm";
import MapView from "./map/MapView";
import LinearView from "./linear/LinearView";
import { Box } from "@mui/joy";
import { STATE_KEY, saveToLocalStorage, exportState } from "./utils/stateUtils";
import { useAppDispatch, useAppSelector } from "@/store/hooks";
import {
  resetSpace,
  setSpace,
  updateView,
  deleteView,
  setDataFiles,
  selectSpace,
  setDependency,
  setConnection,
  deleteConnection,
  deleteDependency,
} from "@/store/features/space/spaceSlice";
import { addLinearGenomeView, addCircosView, addMapView, addCustomMapView } from "./utils/viewUtils";
import FileViewerPanel from "@/components/FileViewerPanel";
import Sidebar from "@/components/Sidebar";
import Multilift from "./multilift/Multilift";
import { getTrackFiles } from "./utils/fileUtils";

interface FileEntry {
  name: string;
  type: string;
  size: number;
}

export default function Page() {

  const dispatch = useAppDispatch();
  const space = useAppSelector(selectSpace);
  const [files, setFiles] = useState<FileEntry[]>([]);
  const [genomeFormOpen, setGenomeFormOpen] = useState(false);

  const fileFormatMapping = {
    "fasta": "karyotype",
    "bed": "line",
    "bedgraph": "line",
    "bedpe": "link",
    "fa": "karyotype",
  }

  const triggerFileRefresh = () => {
    getTrackFiles(space, false).then((trackFiles) => {
      setFiles(trackFiles);
    });
  }

  const [refreshFileViewer, setRefreshFileViewer] = useState(false);

  const deleteTrack = (trackName: string) => {
    setFiles(files.filter((file) => file.name !== trackName));
    dispatch(setDataFiles(space.dataFiles.filter((file) => file !== trackName)));

    const host = process.env.NEXT_PUBLIC_DJANGO_HOST;
    fetch(`${host}/api/timge/delete_track/?uuid=${space.uuid}&track_name=${trackName}`, {
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

  useEffect(() => {
    saveToLocalStorage(space);
  }, [space]);

  useEffect(() => {
    const host = process.env.NEXT_PUBLIC_DJANGO_HOST;
    
    fetch(`${host}/api/timge/get_files_in_path/`, {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({
        uuid: space.uuid,
        path: space.config.working_directory,
      }),
    })
      .then((response) => response.json())
      .then((data) => {
        if (data.status === "success") {
          const trackFiles = data.entries.map((file) => ({
            name: file.name,
            type: file.type,
            size: file.size,
          }));
          setFiles(trackFiles);
        } else {
          console.error("Error fetching tracks:", data.message);
        }
      })
      .catch((err) => {
        console.error("Request failed:", err);
      });
  }, [space.uuid]);

  const uploadTrackFiles = async(files: File[]) => {
    const host = process.env.NEXT_PUBLIC_DJANGO_HOST;
    let formData = new FormData();
    
    files.forEach((track) => {
      formData.append("track_files", track);
    });
    formData.append("uuid", space.uuid);
    fetch(`${host}/api/timge/upload_tracks/`, {
      method: "POST",
      body: formData,
    })
    .then((response) => response.json())
    .then((data) => {
      console.log("Files uploaded successfully", data);
      const trackFiles = files.map((file) => ({
        name: file.name,
        data: file,
        type: file.type,
        size: file.size,
        trackType: fileFormatMapping[file.name.split('.').pop()],
      }));
      getTrackFiles(space, false).then((_files) => {
        setFiles(_files);
      });
    })
    .catch((error) => {
      console.error("Error uploading files", error);
    });
  }

  const addConnection = (viewId: string, linkedViewId: string) => {
    console.log("Adding connection from", viewId, "to", linkedViewId);
    dispatch(setConnection(
      { 
        key: viewId, value: [...(space.connections[viewId] || []), linkedViewId] 
      }));
  }

  const removeConnection = (viewId: string) => {
    const dependents = space.connections[viewId];
    if (dependents) {
      dependents.forEach((dependentId) => {
        dispatch(setDependency(
          { key: dependentId, value: [] }
        ));
      });
    }

    dispatch(deleteConnection(viewId));
  }

  const crossViewActionHandler = (action: string, data: any) => {
    if (action === "generate_heatmap") {
      if (!data || !data.track || !data.reference || !data.segmentA || !data.segmentB || !data.resolution) {
        console.error("Invalid or missing data for generate_heatmap action");
        return;
      }
      console.log("Generating heatmap with data:", data);
      addCustomMapView(dispatch, space,
        {
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
        dispatch(deleteView(viewId));
        dispatch(deleteConnection(viewId));
        dispatch(deleteDependency(viewId));
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
      const dependents = space.connections[viewId];
      if (dependents) {
        dependents.forEach((dependentId) => {
          dispatch(setDependency({ key: dependentId, value: dependencies }));
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
            dispatch(setSpace(state));
          }
        };
        reader.readAsText(file);
      }
    };
    fileInput.click();
  }

  const importTracks = () => {
    setGenomeFormOpen(true);
  }

  const updateViewState = (index: number, updatedConfig: any) => {
    dispatch(updateView({ index, updated: updatedConfig }));
  }

  return <>
    <Header 
      addLinearGenomeView={() => addLinearGenomeView(dispatch, space)}
      addCircosView={() => addCircosView(dispatch, space)}
      addMapView={() => addMapView(dispatch, space)}
      importTracks={importTracks}
      importState={importState}
      exportState={exportState}
      resetState={() => {
        localStorage.removeItem(STATE_KEY);
        setFiles([]);
        dispatch(resetSpace());
      }}
    />
    
    <TrackUploadForm 
      isOpen={genomeFormOpen} 
      onClose={() => setGenomeFormOpen(false)} 
      tracks={files}
      onDeleteTrack={deleteTrack}
      onTrackUpload={(data_files: File[]) => {
        uploadTrackFiles(data_files);
        setRefreshFileViewer(true);
      }}
    />
    { space.config?.multilift_form_open &&
      <Multilift
        triggerFileRefresh={triggerFileRefresh}
      />
    }
    <Box
      sx={{
        display: "flex",
        flexDirection: "row",
        width: "100%",
        height: "100%",
      }}
    >
      {
        space.config?.expanded_sidebar ? (
          <FileViewerPanel
            refreshView={refreshFileViewer}
            setRefreshView={setRefreshFileViewer}
        />
        ) : <Sidebar />
      }

      <Box
        sx={{
          width: space.config?.expanded_sidebar ? "82.5%" : "95%",
          height: "100%",
          backgroundColor: "white",
          padding: "5px",
          overflowY: "auto",
        }}
        >
      <Box
        sx={{
          display: "flex",
          flexDirection: "row",
          flexWrap: "wrap",
          alignItems: "flex-start",
        }}
        >
        {
          space.views.map((view, index) => {
            if (view.type === "linear") {
              return <LinearView 
                        key={index} 
                        trackFiles={files} 
                        viewConfig={view} 
                        handleViewUpdate={updateViewState} 
                        index={index} 
                        crossViewActionHandler={crossViewActionHandler} 
                        dependencies={space.dependencies[view.uuid]}
                        addConnection={addConnection}
                        removeConnection={removeConnection}
                      />  
            }
            else if (view.type === "circos") {
              return <CircosView 
                        key={index}
                        viewConfig={view}
                        handleViewUpdate={updateViewState} 
                        index={index} 
                        crossViewActionHandler={crossViewActionHandler} 
                        dependencies={space.dependencies[view.uuid]}
                        addConnection={addConnection}
                        removeConnection={removeConnection}
                        files={files}
                      />
            }
            else if (view.type === "map") {
              return <MapView 
                        key={index} 
                        trackFiles={files} 
                        viewConfig={view}
                        handleViewUpdate={updateViewState} 
                        index={index} 
                        crossViewActionHandler={crossViewActionHandler} 
                        dependencies={space.dependencies[view.uuid]}
                        addConnection={addConnection}
                        removeConnection={removeConnection}
                      />
            }
          })
        }
        </Box>
      </Box>
    </Box>
    </>;

}