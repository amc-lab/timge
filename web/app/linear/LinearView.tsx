import { useEffect, useState } from "react";
import {
  createViewState,
  JBrowseLinearGenomeView,
} from "@jbrowse/react-linear-genome-view";
import { Box, Select, Option, Typography, Button } from "@mui/joy";
import ParentView from "@/components/ParentView";
import { View } from "../types/state";

interface LinearViewProps {
  trackFiles: any[];
  viewConfig: View;
  handleViewUpdate: (index: number, viewState: View) => void;
  index: number;
  crossViewActionHandler?: any;
  dependencies?: any;
  addConnection?: any;
  removeConnection?: any;
  createdViews: Set<any>;
}

const LinearView = (props: LinearViewProps) => {
  const [reference, setReference] = useState("");
  const [selectedTracks, setSelectedTracks] = useState<string[]>([]);
  const [renderView, setRenderView] = useState(false);
  const [viewState, setViewState] = useState<any>(null);

  const HOST = "http://localhost:3000/";
  const uuid = props.viewConfig.uuid;

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
        name: "Fodor",
        aliases: ["Fodor"],
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
    
    const jbrowseTracks = [
        {
        "trackId": "my_wiggle_track",
          "name": "My Wiggle Track",
          "assemblyNames": ["Fodor"],
          "type": "QuantitativeTrack",
          "adapter": {
            "type": "BedGraphAdapter",
            "bedGraphLocation": {
              "locationType": "UriLocation",
              "uri": "http://localhost:3000/" + selectedTracks[0],
            },
          }
        }
      ];
      
      setViewState(
        createViewState({
          assembly,
          tracks: jbrowseTracks,
        })
      );
    }
  }, [reference, selectedTracks]);

  useEffect(() => {
    if (viewState === null) return;
    console.log("ViewState updated:", viewState);
    const { session } = viewState;
  }, [viewState]);

  return (
    <Box>
      {(!reference || selectedTracks.length === 0 || !renderView) && (
        <ParentView
          viewConfig={props.viewConfig}
          index={props.index}
          crossViewActionHandler={props.crossViewActionHandler}
          >
            <Box sx={{ display: "flex", flexDirection: "column", alignItems: "center" }}>
            <Box sx={{ display: "flex", flexDirection: "row", alignItems: "center" }}>
                <Box
                    sx={{
                        display: "flex",
                        flexDirection: "column",
                        justifyContent: "center",
                        alignItems: "center",
                        marginRight: 2,
                    }}
                >
                    <Typography>Select a reference genome:</Typography>
                    <Select
                    value={reference}
                    onChange={(event, value) => {
                        setReference(value);
                    }}
                    sx={{ width: "200px" }}
                    >
                        {props.trackFiles
                            .filter((file) => file.name.endsWith(".fa") || file.name.endsWith(".fasta"))
                            .map((file, index) => (
                                <Option key={index} value={file.name}>
                                {file.name}
                                </Option>
                        ))}
                    </Select>
                </Box>
                <Box
                    sx={{
                        display: "flex",
                        flexDirection: "column",
                        justifyContent: "center",
                        alignItems: "center",
                        marginRight: 2,
                    }}
                >
                    <Typography>Select tracks to display:</Typography>
                    <Select
                    value={selectedTracks}
                    onChange={(event, value) => {
                        setSelectedTracks(value);
                    }}
                    multiple
                    sx={{ width: "200px"}}
                    >
                        {props.trackFiles
                            .filter((file) => file.name.endsWith(".bedgraph"))
                            .map((file, index) => (
                                <Option key={index} value={file.name}>
                                {file.name}
                                </Option>
                        ))}
                    </Select>
                </Box>
                <Button
                    variant="solid"
                    color="primary"
                    onClick={() => {
                        setRenderView(true);
                    }}
                    >Submit</Button>
                    </Box>
            </Box>

        </ParentView>
      )}

      { renderView && reference && selectedTracks.length > 0 && viewState && (
        <Box sx={{ mt: 4 }}>
          <JBrowseLinearGenomeView viewState={viewState} />
        </Box>
      )}
    </Box>
  );
};

export default LinearView;
