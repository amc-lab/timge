"use client";
import React, { useEffect, useState } from "react";
import {
  Container,
  Box,
  TextField,
  Button,
  Select,
  MenuItem,
  InputLabel,
  FormControl,
  Checkbox,
  ListItemText,
  OutlinedInput,
  Typography,
  Accordion,
  AccordionSummary,
  AccordionDetails,
} from "@mui/material";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";

const Multilift: React.FC = () => {
  const [genomes, setGenomes] = useState({});
  const [genomeInput, setGenomeInput] = useState<string>("");
  const [sequences, setSequences] = useState({});
  const [groups, setGroups] = useState<string[]>(["Group1"]);
  const [liftoverTracks, setLiftoverTracks] = useState<{ file: File; genome: string }[]>([]);

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
    fetch("http://127.0.0.1:8000/multilift/multilift_sequences/", {
      method: "POST",
      body: formData,
    })
      .then((response) => response.json())
      .then((data) => {
        console.log(data);
        for (let i = 0; i < data.length; i++) {
          const seq = data[i];
          const key = [seq[0], seq[1], seq[2]].join(",");

          setSequences((prevSequences) => ({
            ...prevSequences,
            [key]: seq[3],
          }));
        }
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

  const handleLiftoverTrackUpload = (files: FileList | null, genome: string) => {
    if (files) {
      const newTracks = Array.from(files).map((file) => ({ file, genome }));
      setLiftoverTracks((prevTracks) => [...prevTracks, ...newTracks]);
    }
  };

  const generateAlignment = () => {
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

    fetch("http://127.0.0.1:8000/multilift/temp/", {
      method: "POST",
      body: formData,
    })
      .then((response) => {
        if (!response.ok) {
          throw new Error(`HTTP error! status: ${response.status}`);
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
    <Container>
      <Typography variant="h4" gutterBottom>
        Multilift
      </Typography>

      <Accordion defaultExpanded>
        <AccordionSummary expandIcon={<ExpandMoreIcon />}>
          <Typography>Define Genomes</Typography>
        </AccordionSummary>
        <AccordionDetails>
          <Box display="flex" alignItems="center" gap={2}>
            <TextField
              label="Enter genome name"
              value={genomeInput}
              onChange={(e) => setGenomeInput(e.target.value)}
              fullWidth
            />
            <Button variant="contained" onClick={addGenome}>
              Add Genome
            </Button>
          </Box>

          {Object.keys(genomes).length > 0 && (
            <Box mt={2}>
              {Object.keys(genomes).map((genome) => (
                <Box
                  key={genome}
                  display="flex"
                  alignItems="center"
                  justifyContent="space-between"
                  mt={1}
                >
                  <Typography>{genome}</Typography>
                </Box>
              ))}
            </Box>
          )}
        </AccordionDetails>
      </Accordion>

      {Object.keys(genomes).length > 0 && (
        <Accordion defaultExpanded>
          <AccordionSummary expandIcon={<ExpandMoreIcon />}>
            <Typography>Add / Replace Files</Typography>
          </AccordionSummary>
          <AccordionDetails>
            {Object.keys(genomes).map((genome) => (
              <Box key={genome} mt={2}>
                <Typography>{`Upload files for genome: ${genome}`}</Typography>
                <input
                  type="file"
                  onChange={(e) => {
                    if (e.target.files) {
                      addGenomeFile(genome, e.target.files[0]);
                    }
                  }}
                />
              </Box>
            ))}
          </AccordionDetails>
        </Accordion>
      )}

      <Accordion defaultExpanded>
        <AccordionSummary expandIcon={<ExpandMoreIcon />}>
          <Typography>Upload Liftover Tracks</Typography>
        </AccordionSummary>
        <AccordionDetails>
          <Box mt={2}>
            {Object.keys(genomes).map((genome) => (
              <Box key={genome} mt={2}>
                <Typography>{`Upload liftover tracks for genome: ${genome}`}</Typography>
                <input
                  type="file"
                  multiple
                  onChange={(e) => handleLiftoverTrackUpload(e.target.files, genome)}
                />
              </Box>
            ))}
            {liftoverTracks.length > 0 && (
              <Box mt={2}>
                <Typography>Uploaded Liftover Tracks:</Typography>
                <ul>
                  {liftoverTracks.map(({ file, genome }, index) => (
                    <li key={index}>{`${file.name} (Genome: ${genome})`}</li>
                  ))}
                </ul>
              </Box>
            )}
          </Box>
        </AccordionDetails>
      </Accordion>

      {Object.keys(sequences).length > 0 && (
        <Box>
          <Accordion defaultExpanded>
            <AccordionSummary expandIcon={<ExpandMoreIcon />}>
              <Typography>Assign Sequences</Typography>
            </AccordionSummary>
            <AccordionDetails>
              {groups.map((group, index) => (
                <Box key={index} mt={2}>
                  <Typography>{group}</Typography>
                  <FormControl fullWidth>
                    <InputLabel>Select Sequences</InputLabel>
                    <Select
                      multiple
                      value={Object.keys(sequences).filter(
                        (key) => sequences[key] === group
                      )}
                      onChange={(e) => {
                        const selectedKeys = e.target.value as string[];
                        setSequences((prevSequences) => {
                          const updatedSequences = { ...prevSequences };
                          Object.keys(updatedSequences).forEach((key) => {
                            if (selectedKeys.includes(key)) {
                              updatedSequences[key] = group;
                            } else if (updatedSequences[key] === group) {
                              updatedSequences[key] = null;
                            }
                          });
                          return updatedSequences;
                        });
                      }}
                      input={<OutlinedInput label="Select Sequences" />}
                      renderValue={(selected) =>
                        (selected as string[])
                          .map((key) => key.split(",")[2])
                          .join(", ")
                      }
                    >
                      {Object.keys(sequences).map((key) => (
                        <MenuItem key={key} value={key}>
                          <Checkbox checked={sequences[key] === group} />
                          <ListItemText primary={`Seq: ${key.split(",")[2]} | Genome: ${key.split(",")[0]}`} />
                        </MenuItem>
                      ))}
                    </Select>
                  </FormControl>
                </Box>
              ))}

              <Box mt={2} display="flex" justifyContent="space-between">
                <Button variant="contained" onClick={addGroup}>
                  Add Group
                </Button>
                <Button
                  variant="outlined"
                  color="error"
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
                >
                  Remove Group
                </Button>
              </Box>
            </AccordionDetails>
          </Accordion>
          <Box>
            <Button onClick={generateAlignment}>Generate Alignment</Button>
          </Box>
        </Box>
      )}
    </Container>
  );
};

export default Multilift;
