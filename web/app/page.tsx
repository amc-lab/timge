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
import { Track } from "./circos/config/track";
import CircosView from "./circos/CircosView";
import { Box } from "@mui/material";
import { useState } from "react";
import { Assembly } from "./types/genomes";

export default function Page() {
  const viewState = createViewState({
    assembly,
    tracks: tracks_,
  });

  const [tracks, setTracks] = useState<Track[]>([]);
  const [segments, setSegments] = useState<Assembly[]>([]);

  const [components, setComponents] = React.useState([]);
  // const [tracks, setTracks] = useState<Track[]>([]);

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

  return <>
    <Header 
      addLinearGenomeView={addLinearGenomeView}
      addCircosView={addCircosView}
      importTracks={() => {}}
    />
      {
        components.map((component, index) => {
          if (component.viewType === "LinearGenomeView") {
            return <JBrowseLinearGenomeView key={index} viewState={createViewState({
              assembly,
              tracks,
            })} />
          }
          else if (component.viewType === "CircosView") {
            return <CircosView key={index} />
          }
        })
        }
    </>;
}
