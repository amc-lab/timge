"use client";
import React, { useState } from "react";
import { FormControl, Button, Typography, Sheet, Card } from "@mui/joy";
import { TextField } from "@mui/material";
import { Box } from "@mui/system";

const Multilift = () => {
  const [genomes, setGenomes] = useState([{ label: "", file: null }]);

  const handleAddGroup = () => {
    setGenomes([...genomes, { label: "", file: null }]);
  };

  const handleLabelChange = (index, value) => {
    const updatedGenomes = [...genomes];
    updatedGenomes[index].label = value;
    setGenomes(updatedGenomes);
  };

  const handleFileChange = (index, file) => {
    const updatedGenomes = [...genomes];
    updatedGenomes[index].file = file;
    setGenomes(updatedGenomes);
  };

  const handleSubmit = async () => {
    try {
      const response = await fetch("http://127.0.0.1:8000/multilift/", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
      });
  
      if (!response.ok) {
        throw new Error(`Error: ${response.statusText}`);
      }
  
      const result = await response.json();
      console.log("Upload successful:", result);
    } catch (error) {
      console.error("Error during submission:", error);
    }
  };
  
  return (
    <Sheet
      sx={{
        padding: 2,
        minWidth: "100vw",
        minHeight: "100vh",
        display: "flex",
        flexDirection: "column",
        alignItems: "center",
        gap: 2,
      }}
    >
      <Card sx={{ padding: 2, width: "80%" }}>
      <Typography level="h5" sx={{ marginBottom: 2 }}>Sequence Groups</Typography>
        {genomes.map((group, index) => (
          <Box
            key={index}
            sx={{
              display: "flex",
              alignItems: "center",
              gap: 2,
              marginBottom: 2,
            }}
          >
            <TextField
              required
              label={`Sequence Group ${index + 1}`}
              value={group.label}
              onChange={(e) => handleLabelChange(index, e.target.value)}
              sx={{ flex: 1 }}
            />
            <input
              type="file"
              accept=".txt,.fasta,.csv,.fa"
              onChange={(e) =>
                handleFileChange(
                  index,
                  e.target.files ? e.target.files[0] : null
                )
              }
              style={{ flex: 1 }}
            />
          </Box>
        ))}
        <Box sx={{ display: "flex", gap: 2, justifyContent: "center" }}>
            <Button variant="outlined" onClick={handleAddGroup}>
            Add Another Group
            </Button>
            <Button variant="solid" onClick={handleSubmit}>
            Generate Alignment
            </Button>
        </Box>
        </Card>
    </Sheet>
  );
};

export default Multilift;
