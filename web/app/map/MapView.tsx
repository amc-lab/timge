"use client";
import React from "react";
import { View } from "../types/state";
import ParentView from "@/components/ParentView";

interface MapViewProps {
    trackFiles: any[];
    viewConfig: View;
    handleViewUpdate: (index, viewState: View) => void;
    index: number;
}

const MapView = (props: MapViewProps) => {
    return (
        <ParentView
            viewConfig={props.viewConfig}
        >

        </ParentView>
    )
}

export default MapView;