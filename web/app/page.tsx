"use client";
import React, { useEffect } from "react";
import "@fontsource/roboto";
// import {
  // createViewState,
  // JBrowseLinearGenomeView,
// } from "@jbrowse/react-linear-genome-view";
import assembly from "./assembly";
import tracks_ from "./tracks";
import Header from "@/components/Header";
import { v4 as uuidv4 } from "uuid";
import CircosView from "./circos/CircosView";
import { useState } from "react";
import TrackUploadForm from "@/components/TrackUploadForm";
import { State } from "./types/state";
import { ViewType } from "./types/viewTypes";

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
  // const viewState = createViewState({
  //   assembly,
  //   tracks: tracks_,
  // });

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

  const [files, setFiles] = useState<FileEntry[]>([]);
  const [tracks, setTracks] = useState([]);
  const fileInputRef = React.useRef<HTMLInputElement>(null);
  const [genomeFormOpen, setGenomeFormOpen] = useState(false);
  const [spaceState, setSpaceState] = useState<State>(() => loadState());

  useEffect(() => {
    saveState(spaceState);
  }, [spaceState]);

  useEffect(() => {
    const getTrackFiles = () => {
      const host = process.env.NEXT_PUBLIC_DJANGO_HOST;
      fetch(`${host}/api/timge/get_tracks/?uuid=${spaceState.UUID}`, {
        method: "GET",
      })
        .then((response) => response.json())
        .then((data) => {
          if (data.status === "success") {
            const fileMap = data.track_files;
  
            const trackFiles = Object.entries(fileMap).map(([fileName, content]) => ({
              name: fileName,
              data: content,
              trackType: fileFormatMapping[fileName.split('.').pop()],
            }));
  
            setTracks(trackFiles);
            setSpaceState((prevState) => ({
              ...prevState,
              dataFiles: [
                ...prevState.dataFiles,
                ...trackFiles.map((track) => track.name),
              ],
            }));
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

  const addLinearGenomeView = () => {
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
          name: `Linear Genome View ${spaceState.views.length + 1}`,
          description: "Standard linear genome view",
        },
      ],
    });
  }

  const addCircosView = () => {
    const newView = {
      viewType: "CircosView",
    }
    setSpaceState({
      ...spaceState,
      views: [
        ...spaceState.views,
        {
          id: `view-${spaceState.views.length + 1}`,
          type: ViewType.Circos,
          name: `Circos View ${spaceState.views.length + 1}`,
          description: "Circular genome view",
        },
      ],
    });
  }

  const importTracks = () => {
    // fileInputRef.current?.click();
    setGenomeFormOpen(true);
  }

  const handleFileChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    const files = Array.from(event.target.files);
    let circosFiles = [];

    const formData = new FormData();
    formData.append("track_types", JSON.stringify(files.map((file) => fileFormatMapping[file.name.split(".").pop()])));
    files.forEach((file) => {
      formData.append("data_files", new Blob([file]), file.name);
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

  return <>
    <Header 
      addLinearGenomeView={addLinearGenomeView}
      addCircosView={addCircosView}
      importTracks={importTracks}
      importState={() => {
        const state = loadState();
        setSpaceState(state);
      }}
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
      onDeleteTrack={(trackName: string) => {
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
        }
        )
      }}
      onTrackUpload={(files: File[]) => {
        const newTracks = files.map((file) => ({
          name: file.name,
          data: file,
          trackType: fileFormatMapping[file.name.split(".").pop()],
        }));
        setTracks([...tracks, ...newTracks]);
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
    <input type="file" ref={fileInputRef} onChange={handleFileChange} style={{ display: "none" }} />
      {
        spaceState.views.map((view, index) => {
          if (view.type === "linear") {
            return;
            // return <JBrowseLinearGenomeView key={index} viewState={viewState} />
          }
          else if (view.type === "circos") {
            return <CircosView key={index} trackFiles={tracks} />
          }
        })
      }
    </>;
}
