"use client";
import React, { useEffect, useState } from "react";
import {
  Container,
  Box,
  FormControl,
  Typography,
  Snackbar,
  Alert,
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
} from "@mui/joy";
import FormLabel from "@mui/joy/FormLabel";
import Input from "@mui/joy/Input";
import Chip from "@mui/joy/Chip";

const Multilift: React.FC = () => {
  const [alertOpen, setAlertOpen] = useState(false);
  const [alertMessage, setAlertMessage] = useState("Upload successful");
  const [alertSeverity, setAlertSeverity] = useState<
    "error" | "warning" | "info" | "success"
  >("success");
  const [showLoadingBar, setShowLoadingBar] = useState(false);

  const [genomes, setGenomes] = useState({});
  const [genomeInput, setGenomeInput] = useState<string>("");
  const [sequences, setSequences] = useState({});
  const [groups, setGroups] = useState<string[]>(["Group1"]);
  const [liftoverTracks, setLiftoverTracks] = useState<
    { file: File; genome: string }[]
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
        setShowLoadingBar(true);
        const data = await response.json();
        setShowLoadingBar(false);

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
            [key]: seq[3],
          }));
        }
      })
      .catch((error) => {
        console.error("Error fetching sequences:", error);
      });
  };

  useEffect(() => {
    console.log(sequences);
  }, [sequences]);

  const addGroup = () => {
    const newGroup = `Group${groups.length + 1}`;
    setGroups([...groups, newGroup]);
  };

  const assignSequence = (key: string, group: string) => {
    setSequences((prevSequences) => ({
      ...prevSequences,
      [key]: group,
    }));
  };

  const handleLiftoverTrackUpload = (
    files: FileList | null,
    genome: string,
  ) => {
    if (files) {
      const newTracks = Array.from(files).map((file) => ({ file, genome }));
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

    liftoverTracks.forEach(({ file }, index) => {
      formData.append(`uploaded_files`, file);
    });

    const liftoverGenomes = liftoverTracks.map(({ genome }) => genome);
    formData.append("multilift_genomes", JSON.stringify(liftoverGenomes));

    const host = process.env.NEXT_PUBLIC_DJANGO_HOST;
    fetch(`${host}/api/multilift/temp/`, {
      method: "POST",
      body: formData,
    })
      .then((response) => {
        setShowLoadingBar(false);
        if (!response.ok) {
          // setAlertMessage(response.)
          setAlertSeverity("error");
          setAlertOpen(true);
          throw new Error(`HTTP error! status: ${response.status}`);
        } else {
          setAlertMessage("Alignment generated successfully");
          setAlertSeverity("success");
          setAlertOpen(true);
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

  return (
    <Box>
      <LinearProgress
        sx={{
          display: showLoadingBar ? "block" : "none",
          height: "0.5em",
        }}
      />

      <Container
        sx={{
          marginTop: "3em",
        }}
      >
        <Typography variant="h4" gutterBottom>
          Multilift
        </Typography>

        <Box
          sx={{
            marginBottom: "3em",
          }}
        >
          <AccordionGroup>
            <Accordion defaultExpanded>
              <AccordionSummary
                sx={{
                  height: "3em",
                }}
              >
                <Typography>
                  <b>Define Genomes</b>
                </Typography>
              </AccordionSummary>
              <AccordionDetails>
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
                <Box display="flex" flexWrap="wrap" gap={1} mt={1}>
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
                <AccordionSummary>
                  <Typography>
                    <b>Upload Genome</b>
                  </Typography>
                </AccordionSummary>
                <AccordionDetails>
                  {Object.keys(genomes).map((genome) => (
                    <Card
                      key={genome}
                      sx={{
                        padding: "1em",
                        marginBottom: "0.5em",
                      }}
                    >
                      <Typography>{`Upload files for genome: ${genome}`}</Typography>
                      <input
                        type="file"
                        onChange={(e) => {
                          if (e.target.files) {
                            addGenomeFile(genome, e.target.files[0]);
                          }
                        }}
                      />
                    </Card>
                  ))}
                </AccordionDetails>
              </Accordion>
            )}
            {Object.keys(sequences).length > 0 && (
              <Box>
                <Accordion defaultExpanded>
                  <AccordionSummary>
                    <Typography>
                      <b>Upload Data Tracks</b>
                    </Typography>
                  </AccordionSummary>
                  <AccordionDetails>
                    <Box>
                      {Object.keys(genomes).map((genome) => (
                        <Card
                          key={genome}
                          sx={{
                            padding: "1em",
                            marginBottom: "0.5em",
                          }}
                        >
                          <Typography>{`Upload data tracks for genome: ${genome}`}</Typography>
                          <input
                            type="file"
                            multiple
                            onChange={(e) =>
                              handleLiftoverTrackUpload(e.target.files, genome)
                            }
                          />
                        </Card>
                      ))}
                      {liftoverTracks.length > 0 && (
                        <Box mt={2}>
                          <Typography>Uploaded Liftover Tracks:</Typography>
                          <ul>
                            {liftoverTracks.map(({ file, genome }, index) => (
                              <li
                                key={index}
                              >{`${file.name} (Genome: ${genome})`}</li>
                            ))}
                          </ul>
                        </Box>
                      )}
                    </Box>
                  </AccordionDetails>
                </Accordion>
              </Box>
            )}

            {Object.keys(sequences).length > 0 && (
              <Box>
                <Accordion defaultExpanded>
                  <AccordionSummary>
                    <Typography>
                      <b>Assign Sequences</b>
                    </Typography>
                  </AccordionSummary>
                  <AccordionDetails>
                    {groups.map((group, index) => (
                      <Box key={index} mt={2}>
                        <Typography>{group}</Typography>
                        <FormControl fullWidth>
                          {/* <InputLabel>Select Sequences</InputLabel> */}
                          <Select
                            multiple
                            value={Object.keys(sequences).filter(
                              (key) => sequences[key] === group,
                            )}
                            renderValue={(selected) => (
                              <Box sx={{ display: "flex", gap: "0.25rem" }}>
                                {selected.map((selectedOption, index) => (
                                  <Chip
                                    variant="soft"
                                    color="primary"
                                    key={index}
                                  >
                                    {selectedOption.label}
                                  </Chip>
                                ))}
                              </Box>
                            )}
                            onChange={(e, newValue) => {
                              const selectedKeys = newValue as string[];
                              setSequences((prevSequences) => {
                                const updatedSequences = { ...prevSequences };
                                Object.keys(updatedSequences).forEach((key) => {
                                  if (
                                    Array.isArray(selectedKeys) &&
                                    selectedKeys.includes(key)
                                  ) {
                                    updatedSequences[key] = group;
                                  } else if (updatedSequences[key] === group) {
                                    updatedSequences[key] = null;
                                  }
                                });
                                return updatedSequences;
                              });
                            }}
                          >
                            {Object.keys(sequences).map((key) => (
                              <Option key={key} value={key}>
                                {key}
                              </Option>
                            ))}
                          </Select>
                        </FormControl>
                      </Box>
                    ))}

                    <Box mt={2} display="flex" justifyContent="space-between">
                      <Button variant="solid" onClick={addGroup}>
                        Add Group
                      </Button>
                      <Button
                        variant="outlined"
                        onClick={() => {
                          if (groups.length > 1) {
                            const groupToRemove = groups[groups.length - 1];
                            setGroups(groups.slice(0, -1));
                            setSequences((prevSequences) => {
                              const updatedSequences = { ...prevSequences };
                              Object.keys(updatedSequences).forEach((key) => {
                                if (updatedSequences[key] === groupToRemove) {
                                  updatedSequences[key] = null;
                                }
                              });
                              return updatedSequences;
                            });
                          }
                        }}
                        disabled={groups.length <= 1}
                        sx={{ color: "red", borderColor: "red" }}
                      >
                        Remove Group
                      </Button>
                    </Box>
                  </AccordionDetails>
                </Accordion>
                <Box>
                  <Button
                    variant="solid"
                    onClick={generateAlignment}
                    sx={{ width: "100%", height: "3em" }}
                  >
                    Generate Alignment
                  </Button>
                </Box>
              </Box>
            )}
          </AccordionGroup>
        </Box>
        <Snackbar
          key={"right"}
          open={alertOpen}
          anchorOrigin={{ horizontal: "right", vertical: "bottom" }}
          autoHideDuration={6000}
        >
          <Alert
            severity={alertSeverity}
            sx={{
              width: "100%",
              fontSize: "1em",
            }}
          >
            {alertMessage}
          </Alert>
        </Snackbar>
      </Container>
    </Box>
  );
};

export default Multilift;
