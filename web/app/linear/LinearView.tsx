"use client";
import { useEffect, useRef, useState, StrictMode } from "react";
import {
  createViewState,
  JBrowseLinearGenomeView,
} from "@jbrowse/react-linear-genome-view";
import { Box, Select, Option, Typography, Button } from "@mui/joy";
import ParentView from "@/components/ParentView";
import { View } from "@/store/features/views/types";
import { useAppSelector, useAppDispatch} from "@/store/hooks";
import TrackSelector from "./components/TrackSelector";
import IGVBrowser from "./IGVBrowser";
import { updateViewProps, updateViewConfig} from "@/store/features/space/spaceSlice";

interface LinearViewProps {
  trackFiles: any[];
  viewConfig: View;
  handleViewUpdate: (index: number, viewState: View) => void;
  index: number;
  crossViewActionHandler?: any;
  dependencies?: any;
  addConnection?: any;
  removeConnection?: any;
  // createdViews: Set<any>;
}

const LinearView = (props: LinearViewProps) => {
  const space = useAppSelector((state) => state.space);
  const [reference, setReference] = useState(
    space.views.find(v => v.uuid === props.viewConfig.uuid).config.reference || ""
  );
  const [selectedTracks, setSelectedTracks] = useState<string[]>(
    space.views.find(v => v.uuid === props.viewConfig.uuid).config.trackFiles || []
  );
  const [renderView, setRenderView] = useState(false);
  const [viewState, setViewState] = useState<any>(null);
  const [openTrackSelector, setOpenTrackSelector] = useState(false);

  const dispatch = useAppDispatch();
  

  const HOST = "https://timge.doc.ic.ac.uk/uploads/" + space.uuid + "/";
  
  useEffect(() => {
    console.log("ViewState updated:", viewState);
    // const { displayedRegions, bpPerPx, offsetPx } = viewState;
    // console.log("displayedRegions", displayedRegions);
    // console.log("bpPerPx", bpPerPx);
    // console.log("offsetPx", offsetPx);
  }
  , [viewState]);
  
//   // Filter reference tracks (e.g., .fasta/.fa)
//   const referenceTracks = props.trackFiles.filter((file) =>
//     file.endsWith(".fa") || file.endsWith(".fasta")
//   );

//   // Filter data tracks (e.g., .bedgraph/.bigwig/.gff etc.)
//   const dataTracks = props.trackFiles.filter(
//     (file) =>
//       !file.endsWith(".fa") &&
//       !file.endsWith(".fasta")
//   );

  useEffect(() => {
    if (reference && selectedTracks.length > 0) {
    //   const assembly = {
    //     name: "customAssembly",
    //     sequence: {
    //       type: "ReferenceSequenceTrack",
    //       trackId: "ref-track",
    //       adapter: {
    //         type: "IndexedFastaAdapter",
    //         fastaLocation: {
    //           uri: reference,
    //           locationType: "UriLocation",
    //         },
    //         faiLocation: {
    //           uri: `${reference}.fai`,
    //           locationType: "UriLocation",
    //         },
    //       },
    //     },
    //   };

    //   const jbrowseTracks = selectedTracks.map((track: any, i: number) => ({
    //     type: "QuantitativeTrack",
    //     trackId: `track-${i}`,
    //     name: track.name || `Track ${i + 1}`,
    //     assemblyNames: ["customAssembly"],
    //     adapter: {
    //       type: "BedGraphAdapter",
    //       bedGraphLocation: {
    //         uri: track.uri,
    //         locationType: "UriLocation",
    //       },
    //       index: {
    //         location: {
    //           uri: `${track.uri}.tbi`,
    //           locationType: "UriLocation",
    //         },
    //       },
    //     },
    //   }));
    const assembly = {
        name: reference,
        aliases: [reference],
        sequence: {
          type: "ReferenceSequenceTrack",
          trackId: "GRCh38-ReferenceSequenceTrack",
          adapter: {
            type: "IndexedFastaAdapter",
            fastaLocation: {
              uri: HOST + reference,
            },
            faiLocation: {
              uri: HOST + reference + ".fai",
            }
          },
        },
      };
    
    // const jbrowseTracks = [
    //     {
    //     "trackId": selectedTracks[0],
    //       "name": selectedTracks[0],
    //       "assemblyNames": [reference],
    //       "type": "QuantitativeTrack",
    //       "adapter": {
    //         "type": "BedGraphAdapter",
    //         "bedGraphLocation": {
    //           "locationType": "UriLocation",
    //           "uri": HOST + selectedTracks[0],
    //         },
    //       }
    //     }
    //   ];

    const jbrowseTracks = selectedTracks.map((track) => ({
      trackId: track,
      name: track,
      assemblyNames: [reference],
      type: "QuantitativeTrack",
      adapter: {
      type: "BedGraphAdapter",
      bedGraphLocation: {
        locationType: "UriLocation",
        uri: HOST + track,
      },
      },
    }));
      
      setViewState(
        createViewState({
          assembly,
          tracks: jbrowseTracks,
        })
      );
    }
  }, [reference, selectedTracks]);

  const generateFai = async() => {
    const host = process.env.NEXT_PUBLIC_DJANGO_HOST;
    const queryParams = new URLSearchParams({
      uuid: props.viewConfig.uuid,
      genome_path: reference,
    }).toString();
    await fetch(`${host}/api/timge/generate_fai/?${queryParams}`, {
      method: "GET",
      headers: {
        "Content-Type": "application/json",
      },
    })
      .then((response) => response.json())
      .then((data) => {
        if (data.status === "success") {
          console.log("FAI file generated successfully");
          return true;
        } else {
          console.error("Failed to fetch segments", data.message);
          return false;
        }
      });
  }

  useEffect(() => {
    if (viewState === null) return;
    console.log("ViewState updated:", viewState);
    const { session } = viewState;
  }, [viewState]);

  return (
    <Box sx={{
      display: "flex",
      width: "100%",
    }}
    >
      <ParentView
          viewConfig={props.viewConfig}
          index={props.index}
          userActions={{
            "Clear Tracks": () => {
              dispatch(updateViewConfig({
                uuid: props.viewConfig.uuid,
                config: {
                    ...props.viewConfig.config,
                    reference: "",
                    trackFiles: []
                }}))
              setReference("");
              setSelectedTracks([]);
          },
            "Select Tracks": () => {
              setOpenTrackSelector(true);
            }
          }}
        >
      {(!reference) && (
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
              variant="solid"
              color="primary"
              onClick={() => setOpenTrackSelector(true)}
              sx={{
                marginTop: "10px",
              }}
            >
              Select Tracks
            </Button>
        </Box>
      )}

      {openTrackSelector && (
        <TrackSelector
          onClose={() => setOpenTrackSelector(false)}
          onConfirm={(referencePath, selectedTrackPaths) => {
            setReference(referencePath);
            setSelectedTracks(selectedTrackPaths);
            setOpenTrackSelector(false);
            setRenderView(true);
            dispatch(updateViewConfig({
              uuid: props.viewConfig.uuid,
              config: {
                  ...props.viewConfig.config,
                  reference: referencePath,
                  trackFiles: selectedTrackPaths
              }}))
          }}
        />
      )}

      { reference && (
        <Box sx={{ width: "100%" }}>
          {/* <JBrowseLinearGenomeView 
          viewState={viewState} 
          /> */}
          <IGVBrowser
            reference={reference}
            trackFiles={selectedTracks}
          />
        </Box>
      )}
    </ParentView>
    </Box>
  );
};

export default LinearView;
