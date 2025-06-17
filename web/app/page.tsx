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
  setDiffStructureFormOpen,
} from "@/store/features/space/spaceSlice";
import { addLinearGenomeView, addCircosView, addMapView, addCustomMapView } from "./utils/viewUtils";
import FileViewerPanel from "@/components/FileViewerPanel";
import Sidebar from "@/components/Sidebar";
import Multilift from "./multilift/Multilift";
import { getTrackFiles } from "./utils/fileUtils";
import { fetchFiles } from "@/store/features/files/fileSlice";
import { File as FileType } from "@/store/features/files/types";
import LinearProgress from '@mui/material/LinearProgress';
import { setLoading } from "@/store/features/site/siteSlice";
import ShapeReactivityForm from "./diffStructure/Form";
import { store } from "@/store";

interface FileEntry {
  name: string;
  type: string;
  size: number;
}

export default function Page() {

  const dispatch = useAppDispatch();
  const space = useAppSelector(selectSpace);
  // const [files, setFiles] = useState<FileEntry[]>([]);
  const [files, setFiles] = useState<FileType[]>([]);
  const [genomeFormOpen, setGenomeFormOpen] = useState(false);
  const _files = useAppSelector((state) => state.files);
  const isLoading = useAppSelector((state) => state.site.isLoading);

  useEffect(() => {
    if (space.uuid) {
      const fetchData = async () => {
        const _files = await dispatch(fetchFiles({ uuid: space.uuid, path: [] }));
        return _files;
      };
      fetchData().then((result) => {
        if (result.meta.requestStatus === "fulfilled") {
          const files = result.payload as FileType[];
          setFiles(files);
        }
      });
    }
  }, [space.uuid]);
  
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
        dispatch(fetchFiles({ uuid: space.uuid, path: [] }));
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
    // fetch(`${host}/api/timge/upload_tracks/`, {
    //   method: "POST",
    //   body: formData,
    // })
    // .then((response) => response.json())
    // .then((data) => {
    //   console.log("Files uploaded successfully", data);
    //   dispatch(fetchFiles({ uuid: space.uuid, path: [] }));
    // })
    // .catch((error) => {
    //   console.error("Error uploading files", error);
    // });

    try {
      const response = await fetch(`${host}/api/timge/upload_tracks/`, {
        method: "POST",
        body: formData,
      });
  
      const data = await response.json();
      console.log("Files uploaded successfully", data);
  
      await dispatch(fetchFiles({ uuid: space.uuid, path: [] }));
    } catch (error) {
      console.error("Error uploading files", error);
    } finally {
      dispatch(setLoading(false));
    }
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
    else if (action === "propagate_loci") {
      console.log("Propagating loci update:", data);
      const currentSpace = store.getState().space;
      const { viewId, loci } = data;
      console.log("View ID:", viewId, "Loci:", loci);
      console.log("Space connections:", currentSpace.connections);
      const dependents = currentSpace.connections[viewId];
      console.log(dependents);
      if (dependents) {
        dependents.forEach((dependentId) => {
          dispatch(setDependency({ key: dependentId, value: loci }));
        });
      }
    }
    else if (action === "update_view") {
      const { index, updated } = data;
      updateViewState(index, updated);
    }
    else if (action === "set_diff_structure_form_open") {
      dispatch(setDiffStructureFormOpen(data.open));
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
  {isLoading &&
  <LinearProgress
    sx={{
      width: "100%",
      position: "absolute",
      top: 0,
      left: 0,
      right: 0,
      zIndex: 2000,
    }}
    ></LinearProgress>
  }
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
        ) : space.config?.diff_structure_form_open ? (
          <ShapeReactivityForm
            open={space.config?.diff_structure_form_open}
            onClose={() => dispatch(setDiffStructureFormOpen(false))}
            onSubmit={(conditions) => {
              console.log("Conditions submitted:", conditions);
              dispatch(setSpace({ ...space, config: { ...space.config, diff_structure_form_open: false } }));
            }}
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
          // alignItems: "flex-end"
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