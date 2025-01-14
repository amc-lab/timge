"use client";
import React, { useState } from "react";
import { FormControl, Button, Typography, Sheet, Card } from "@mui/joy";
import { TextField, Select, MenuItem } from "@mui/material";
import { Box } from "@mui/system";

const Multilift = () => {
  const [genomes, setGenomes] = useState([{ name: "", referenceTrack: null }]);
  const [aligner, setAligner] = useState("mafft");
  const [outputFormat, setOutputFormat] = useState(".zip");

  const handleAddGenome = () => {
    setGenomes([...genomes, { name: "", referenceTrack: null }]);
  };

  const handleNameChange = (index, value) => {
    const updatedGenomes = [...genomes];
    updatedGenomes[index].name = value;
    setGenomes(updatedGenomes);
  };

  const handleFileChange = (index, file) => {
    const updatedGenomes = [...genomes];
    updatedGenomes[index].referenceTrack = file;
    setGenomes(updatedGenomes);
  };

  const handleSubmit = async () => {
    const formData = new FormData();

    formData.append("genomes", JSON.stringify(genomes.map((genome) => genome.name)));

    genomes.forEach((genome) => {
      if (genome.referenceTrack) {
        formData.append("sequences", genome.referenceTrack);
      }
    });

    formData.append("aligner", JSON.stringify(aligner));
    formData.append("output_format", JSON.stringify(outputFormat));

    try {
      const response = await fetch("http://127.0.0.1:8000/multilift/generate_alignment/", {
        method: "POST",
        body: formData,
      });

      if (!response.ok) {
        throw new Error(`Error: ${response.statusText}`);
      }

      const blob = await response.blob();
      const url = window.URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.style.display = "none";
      a.href = url;
      if (outputFormat === ".tar.gz")
        a.download = "multilift.tar.gz";
      else
        a.download = "multilift.zip";
      document.body.appendChild(a);
      a.click();
      window.URL.revokeObjectURL(url);
    } catch (error) {
      console.error("Error during submission:", error);
    }
  };

  return (
    <Sheet
      sx={{
        padding: 4,
        minWidth: "100vw",
        minHeight: "100vh",
        display: "flex",
        flexDirection: "column",
        alignItems: "center",
        gap: 3,
      }}
    >
      <Card sx={{ padding: 3, width: "90%", maxWidth: "600px" }}>
        {genomes.map((genome, index) => (
          <Box
            key={index}
            sx={{
              display: "flex",
              flexDirection: "column",
              gap: 2,
              marginBottom: 3,
            }}
          >
            <TextField
              required
              label={`Genome Name ${index + 1}`}
              value={genome.name}
              onChange={(e) => handleNameChange(index, e.target.value)}
              sx={{ width: "100%" }}
            />
            <input
              type="file"
              accept=".txt,.fasta,.csv,.fa"
              onChange={(e) =>
                handleFileChange(index, e.target.files ? e.target.files[0] : null)
              }
            />
          </Box>
        ))}
        <Box sx={{ display: "flex", flexDirection: "column", gap: 2 }}>
          <FormControl>
            <Typography>Aligner</Typography>
            <Select
              value={aligner}
              onChange={(e) => setAligner(e.target.value)}
            >
              <MenuItem value="mafft">MAFFT</MenuItem>
              <MenuItem value="clustalo">Clustal</MenuItem>
              <MenuItem value="muscle">MUSCLE</MenuItem>
              <MenuItem value="kalign">Kalign</MenuItem>
            </Select>
          </FormControl>
          <FormControl>
            <Typography>Output Format</Typography>
            <Select
              value={outputFormat}
              onChange={(e) => setOutputFormat(e.target.value)}
            >
              <MenuItem value=".zip">zip</MenuItem>
              <MenuItem value=".tar.gz">tar.gz</MenuItem>
            </Select>
          </FormControl>
        </Box>
        <Box sx={{ display: "flex", gap: 2, justifyContent: "center", marginTop: 3 }}>
          <Button variant="outlined" onClick={handleAddGenome}>
            Add Genome
          </Button>
          <Button variant="solid" onClick={handleSubmit}>
            Create Alignment
          </Button>
        </Box>
      </Card>
    </Sheet>
  );
};

export default Multilift;
