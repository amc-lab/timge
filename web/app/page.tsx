"use client";
import React from "react";
import "@fontsource/roboto";
import {
  createViewState,
  JBrowseLinearGenomeView,
} from "@jbrowse/react-linear-genome-view";
import assembly from "./assembly";
import tracks_ from "./tracks";
import Button from "@mui/material/Button";
import Header from "@/components/Header";
import { Track, TrackType } from "./circos/config/track";
import CircosView from "./circos/CircosView";
import { Box } from "@mui/material";
import { useState } from "react";
import { Assembly } from "./types/genomes";

interface FileEntry {
  name: string;
  data: any;
  trackType: string;
}
export default function Page() {
  const viewState = createViewState({
    assembly,
    tracks: tracks_,
  });

  const fileFormatMapping = {
    "fasta": "karyotype",
    "bed": "line",
    "bedgraph": "line",
    "bedpe": "link",
    "fa": "karyotype",
  }

  const [files, setFiles] = useState<FileEntry[]>([]);
  const [tracks, setTracks] = useState([]);
  const [components, setComponents] = React.useState([]);
  const fileInputRef = React.useRef<HTMLInputElement>(null);

  const addLinearGenomeView = () => {
    const newView = {
      viewType: "LinearGenomeView",
    }
    setComponents([...components, newView]);
  }

  const addCircosView = () => {
    const newView = {
      viewType: "CircosView",
    }
    setComponents([...components, newView]);
  }

  const importTracks = () => {
    fileInputRef.current?.click();
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
    />
    <input type="file" ref={fileInputRef} onChange={handleFileChange} style={{ display: "none" }} />
      {
        components.map((component, index) => {
          if (component.viewType === "LinearGenomeView") {
            return <JBrowseLinearGenomeView key={index} viewState={viewState} />
          }
          else if (component.viewType === "CircosView") {
            return <CircosView key={index} trackFiles={tracks} />
          }
        })
        }
    </>;
}
