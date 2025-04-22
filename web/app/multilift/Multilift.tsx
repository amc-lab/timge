"use client";
import React, { useEffect, useState, useCallback } from "react";
import {
  Container,
  Box,
  FormControl,
  Typography,
  Snackbar,
  Alert,
  Button as MuiButton,
  LinearProgress,
} from "@mui/material";
import Button from "@mui/joy/Button";
import Card from "@mui/joy/Card";
import {
  Accordion,
  AccordionSummary,
  AccordionDetails,
  AccordionGroup,
  Option,
  Select,
  IconButton,
} from "@mui/joy";
import FormLabel from "@mui/joy/FormLabel";
import Input from "@mui/joy/Input";
import Chip from "@mui/joy/Chip";
import GenomeFileUploadBox from "./components/GenomeFileUpload";
import DataTrackFileUploadBox from "./components/DataTrackFileUpload";
import ChipDelete from '@mui/joy/ChipDelete';
import DataTrackSelect from "./components/DataTrackSelect";
import CloseIcon from '@mui/icons-material/Close';
import { useAppDispatch, useAppSelector } from "@/store/hooks";
import { setMultiliftFormOpen } from '@/store/features/space/spaceSlice';

interface MultiliftProps {
  triggerFileRefresh: () => void;
}

const Multilift: React.FC<MultiliftProps> = ({triggerFileRefresh}) => {
  const [showLoadingBar, setShowLoadingBar] = useState(false);

  const [genomes, setGenomes] = useState({});
  const [genomeInput, setGenomeInput] = useState<string>("");
  const [sequences, setSequences] = useState({});
  const [groups, setGroups] = useState<string[]>(["Consensus"]);
  const [liftoverTracks, setLiftoverTracks] = useState<
    { file: File; genome: string | null }[]
  >([]);

  const addGenome = () => {
    setGenomes({
      ...genomes,
      [genomeInput]: null,
    });
    setGenomeInput("");
  };

  const addGenomeFile = (genome: string, file: File) => {
    setGenomes({
      ...genomes,
      [genome]: file,
    });

    const formData = new FormData();
    formData.append("genome", genome);
    formData.append("genome_file", file);

    const host = process.env.NEXT_PUBLIC_DJANGO_HOST;

    fetch(`${host}/api/multilift/multilift_sequences/`, {
      method: "POST",
      body: formData,
    })
      .then(async (response) => {
        const data = await response.json();

        if (data.message) {
        }
        console.log("Response Message:", data.message);

        if (!Array.isArray(data)) {
          console.error("Unexpected response format:", data);
          return;
        }

        for (let i = 0; i < data.length; i++) {
          const seq = data[i];
          const key = [seq[0], seq[1], seq[2]].join(",");

          setSequences((prevSequences) => ({
            ...prevSequences,
            [key]: "Consensus",
          }));
        }
      })
      .catch((error) => {
        console.error("Error fetching sequences:", error);
      });
  };

  const removeGenomeFile = (genome: string) => () => {
    setGenomes({
      ...genomes,
      [genome]: null,
    });
  };

  const removeLiftoverTrack = (index: number) => {
    setLiftoverTracks((prevTracks) => {
      const newTracks = [...prevTracks];
      newTracks.splice(index, 1);
      return newTracks;
    });
  }

  const handleLiftoverTrackSelection = (file: File, genome: string) => {
    setLiftoverTracks((prevTracks) => {
      const newTracks = [...prevTracks];
      const index = newTracks.findIndex((track) => track.file.name === file.name);
      if (index >= 0) {
        newTracks[index] = { file, genome };
      }
      return newTracks;
    });
  };

  useEffect(() => {
    console.log(sequences);
  }, [sequences]);

  const handleLiftoverTrackUpload = (
    files: FileList | null,
  ) => {
    if (files) {
      const newTracks = Array.from(files).map((file) => ({ file, genome: null }));
      setLiftoverTracks((prevTracks) => [...prevTracks, ...newTracks]);
    }
  };

  const generateAlignment = () => {
    setShowLoadingBar(true);
    const formData = new FormData();

    formData.append("genomes", JSON.stringify(Object.keys(genomes)));

    Object.entries(genomes).forEach(([genome, file]) => {
      if (file instanceof File) {
        formData.append("genome_files", file);
      }
    });
    formData.append("sequences", JSON.stringify(sequences));
    formData.append("groups", JSON.stringify(groups));
    formData.append("aligner", "mafft");
    formData.append("download_format", ".zip");

    liftoverTracks.forEach(({ file, genome }, index) => {
      formData.append(`uploaded_files`, file);
    });

    const uploadGenomes = [];
    for (let i = 0; i < liftoverTracks.length; i++) {
      const { file, genome } = liftoverTracks[i];
      if (genome) {
        formData.append(`liftover_files`, file);
        uploadGenomes.push(genome);
      }
    }
    formData.append("multilift_genomes", JSON.stringify(uploadGenomes));

    const host = process.env.NEXT_PUBLIC_DJANGO_HOST;
    fetch(`${host}/api/multilift/temp/`, {
      method: "POST",
      body: formData,
    })
      .then((response) => {
        setShowLoadingBar(false);
        if (!response.ok) {
          throw new Error(`HTTP error! status: ${response.status}`);
        } else {
        }
        return response.blob();
      })
      .then((blob) => {
        const url = window.URL.createObjectURL(blob);
        const a = document.createElement("a");
        a.href = url;
        a.download = "multilift.zip";
        document.body.appendChild(a);
        a.click();
        a.remove();
        window.URL.revokeObjectURL(url);
      })
      .catch((error) => {
        console.error("Error generating alignment:", error);
      });
  };

  const space = useAppSelector((state) => state.space);
  const dispatch = useAppDispatch();

  const saveAlignment = () => {
    const formData = new FormData();

    formData.append("genomes", JSON.stringify(Object.keys(genomes)));

    Object.entries(genomes).forEach(([genome, file]) => {
      if (file instanceof File) {
        formData.append("genome_files", file);
      }
    });
    formData.append("sequences", JSON.stringify(sequences));
    formData.append("groups", JSON.stringify(groups));
    formData.append("aligner", "mafft");
    formData.append("download_format", ".zip");
    formData.append("uuid", space.uuid);

    liftoverTracks.forEach(({ file, genome }, index) => {
      formData.append(`uploaded_files`, file);
    });

    const uploadGenomes = [];
    for (let i = 0; i < liftoverTracks.length; i++) {
      const { file, genome } = liftoverTracks[i];
      if (genome) {
        formData.append(`liftover_files`, file);
        uploadGenomes.push(genome);
      }
    }
    formData.append("multilift_genomes", JSON.stringify(uploadGenomes));

    const host = process.env.NEXT_PUBLIC_DJANGO_HOST;
    fetch(`${host}/api/multilift/temp/`, {
      method: "POST",
      body: formData,
    })
      .then((response) => {
        setShowLoadingBar(false);
        if (!response.ok) {
          throw new Error(`HTTP error! status: ${response.status}`);
        }
        return response.blob();
      })
      .catch((error) => {
        console.error("Error generating alignment:", error);
      });
    
      triggerFileRefresh();
  }


  return (
    <Box
    sx={{
        position: "fixed",
        top: 0,
        left: 0,
        width: "100vw",
        height: "100vh",
        backgroundColor: "rgba(0, 0, 0, 0.5)", // semi-transparent black
        zIndex: 1299, // high z-index to overlay everything
        display: "flex",
        justifyContent: "center",
        alignItems: "center",
    }}
    >
      <Box
          sx={{
          width: "80vw",
          maxHeight: "90vh",
          backgroundColor: "white",
          borderRadius: "1em",
          overflowY: "auto",
          boxShadow: 24,
          padding: "2em",
          }}
      >
        <LinearProgress
        sx={{
            display: showLoadingBar ? "block" : "none",
            height: "0.5em",
        }}
        />

        <Box
          sx={{
            width: "100%",
            display: "flex",
            justifyContent: "space-between",
            alignItems: "center",
          }}
          >
          <Typography variant="h4" gutterBottom>
          Multilift
          </Typography>
          <IconButton
              variant="plain"
              size="sm"
              onClick={() => dispatch(setMultiliftFormOpen(false))}
          >
              <CloseIcon />
          </IconButton>
        </Box>

        <Box sx={{ marginBottom: "1em" }}>
          <AccordionGroup>
            <Accordion defaultExpanded>
              <AccordionSummary
                sx={{
                  height: "4em",
                }}
              >
                <Typography>
                  <b>Provide Genome Names</b>
                </Typography>
              </AccordionSummary>
              <AccordionDetails
              >
                <Box display="flex" alignItems="center" gap={2}>
                    <Box
                    display="block"
                    sx={{
                      width: "85%",
                      paddingBottom: "1em",
                    }}
                    >
                    <FormLabel>Genome</FormLabel>
                    <Input
                      placeholder="Enter here..."
                      value={genomeInput}
                      onChange={(e) => setGenomeInput(e.target.value)}
                      onKeyDown={(e) => {
                      if (e.key === "Enter") {
                        addGenome();
                      }
                      }}
                      sx={{
                      height: "3em",
                      }}
                    />
                    </Box>

                  <Button
                    variant="solid"
                    onClick={addGenome}
                    sx={{
                      width: "15%",
                      height: "3.25em",
                    }}
                  >
                    Add Genome
                  </Button>
                </Box>

                <Typography>Genomes Added:</Typography>
                <Box display="flex" flexWrap="wrap" gap={1} mt={1}
                  sx={{
                    paddingBottom: "1em",
                  }}
                >
                  {Object.keys(genomes).map((genome, index) => (
                    <Chip key={index} variant="soft" color="primary">
                      {genome}
                    </Chip>
                  ))}
                </Box>
              </AccordionDetails>
            </Accordion>

            {Object.keys(genomes).length > 0 && (
              <Accordion defaultExpanded>
                <AccordionSummary
                  sx={{
                    height: "4em",
                  }}
                >
                  <Typography>
                    <b>Upload Genome</b>
                  </Typography>
                </AccordionSummary>
                <AccordionDetails>
                  <Box
                    sx={{
                      display: "flex",
                      flexWrap: "wrap",
                      gap: "1em",
                    }}
                  >
                  {Object.keys(genomes).map((genome) => (
                    <GenomeFileUploadBox
                      key={genome}
                      genome={genome}
                      onGenomeFileUpload={(file) => addGenomeFile(genome, file)}
                      onGenomeFileRemove={removeGenomeFile(genome)}
                    />
                  ))}
                  </Box>
                </AccordionDetails>
              </Accordion>
            )}
            {Object.keys(sequences).length > 0 && (
              <Box>
                <Accordion defaultExpanded>
                  <AccordionSummary
                    sx={{
                      height: "4em",
                    }}
                  >
                    <Typography>
                      <b>Upload Data Tracks</b>
                    </Typography>
                  </AccordionSummary>
                  <AccordionDetails>
                    <Box>
                      <Box
                      sx={{
                        display: "flex",
                        flexWrap: "wrap",
                        gap: "1em",
                      }}
                    >
                        <>
                          <DataTrackFileUploadBox
                            onDataTrackFileUpload={(fileList) =>
                              handleLiftoverTrackUpload(fileList)
                            }
                          />
                        </>
                    </Box>
                      {liftoverTracks.length > 0 && (
                        <Box mt={2}
                          sx={{
                            paddingBottom: "1em",
                          }}
                        >
                          <Typography>Uploaded Liftover Tracks:</Typography>
                          <ul>
                            {liftoverTracks.map(({ file, genome }, index) => (
                              <Chip
                              key={index}
                              variant="soft"
                              color="primary"
                              sx={{ margin: "0.25em" }}
                              endDecorator={<ChipDelete onDelete={() => removeLiftoverTrack(index)} />}
                              >
                              {file.name}
                              </Chip>
                            ))}
                          </ul>
                        </Box>
                      )}
                    </Box>
                    <Box sx={{                       
                      display: "flex",
                      flexWrap: "wrap",
                      gap: "1em", 
                    }}
                    >
                      {
                        Object.keys(genomes)
                          .map((genome, index) => (
                            <DataTrackSelect
                              key={index}
                              dataTracks={liftoverTracks
                                .filter((track) => track.genome === null || track.genome === genome)
                                .map((track) => track.file)}
                              genome={genome}
                              setSelectedDataTrack={handleLiftoverTrackSelection}
                            />
                        ))
                      }
                    </Box>
                  </AccordionDetails>
                </Accordion>
              </Box>
            )}

            {Object.keys(sequences).length > 0 && (
              <Box>
                <Box
                  sx={{
                    display: "flex",
                    alignContent: "center",
                    justifyContent: "center",
                    marginTop: "1em",
                  }}  
                >
                  <Button
                    variant="solid"
                    onClick={generateAlignment}
                    sx={{ minWidth: "25%", height: "3em", marginRight: "1em" }}
                  >
                    Download Alignment
                  </Button>
                  <Button
                    variant="solid"
                    onClick={saveAlignment}
                    sx={{ minWidth: "25%", height: "3em", marginLeft: "1em" }}
                  >
                    Save Alignment to Workspace
                  </Button>
                </Box>
              </Box>
            )}
          </AccordionGroup>
        </Box>
      </Box>
    </Box>
  );
};

export default Multilift;
