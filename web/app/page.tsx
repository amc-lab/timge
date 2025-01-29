"use client";
import React from "react";
import "@fontsource/roboto";
import {
  createViewState,
  JBrowseLinearGenomeView,
} from "@jbrowse/react-linear-genome-view";
import assembly from "./assembly";
import tracks from "./tracks";

export default function Page() {
  const viewState = createViewState({
    assembly,
    tracks,
  });

  return <>{/* <JBrowseLinearGenomeView viewState={viewState} /> */}</>;
}
