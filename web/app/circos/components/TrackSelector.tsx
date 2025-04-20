import React, { useState } from "react";
import {
  Box,
  Checkbox,
  Typography,
  IconButton,
  Sheet,
  Button,
  Divider,
  Card,
} from "@mui/joy";
import ArrowUpwardIcon from "@mui/icons-material/ArrowUpward";
import ArrowDownwardIcon from "@mui/icons-material/ArrowDownward";
import { Track } from "../config/track";

interface TrackSelectorProps {
  // tracks: Track[];
  // trackFiles: any[];
  files: any[];
  onClose: () => void;
  onConfirm: (selectedTracks: Track[]) => void;
}

const TrackSelector: React.FC<TrackSelectorProps> = ({
  // tracks,
  // trackFiles,
  files,
  onClose,
  onConfirm,
}) => {
  const [selectedTracks, setSelectedTracks] = useState<Track[]>([]);

  const toggleTrack = (track: Track) => {
    setSelectedTracks((prev) => {
      if (prev.includes(track)) {
        return prev.filter((t) => t !== track);
      } else {
        return [...prev, track];
      }
    });
  };

  const handleMove = (index: number, direction: number) => {
    const newIndex = index + direction;
    if (newIndex < 0 || newIndex >= selectedTracks.length) return;
    const newOrder = [...selectedTracks];
    const [movedItem] = newOrder.splice(index, 1);
    newOrder.splice(newIndex, 0, movedItem);
    setSelectedTracks(newOrder);
  };

  const handleConfirm = () => {
    onConfirm(selectedTracks);
  };

  return (
    <Box
      sx={{
        position: "fixed",
        top: 0,
        left: 0,
        width: "100vw",
        height: "100vh",
        backgroundColor: "rgba(0, 0, 0, 0.6)",
        display: "flex",
        justifyContent: "center",
        alignItems: "center",
        zIndex: 1300,
      }}
    >
      <Sheet
        sx={{
          width: "80vw",
          maxHeight: "80vh",
          backgroundColor: "white",
          borderRadius: "12px",
          padding: 4,
          overflowY: "auto",
        }}
      >
        <Typography level="h4" mb={2}>
          Track Selector
        </Typography>
        <Divider sx={{ mb: 2 }} />

        {/* Box 1: All tracks with checkboxes */}
        <Typography mb={1}>
          Available Tracks
        </Typography>
        <Box
          sx={{
            display: "flex",
            flexDirection: "column",
            gap: 1,
            mb: 3,
            p: 2,
            border: "1px solid #ccc",
            borderRadius: "8px",
            backgroundColor: "#f9f9f9",
          }}
        >
          {files.map((track, index) => (
            <Box
              key={index}
              sx={{
                display: "flex",
                alignItems: "center",
                gap: 1,
              }}
            >
              <Checkbox
                checked={selectedTracks.includes(track)}
                onChange={() => toggleTrack(track)}
              />
              <Typography>{track.name}</Typography>
              {/* <Typography>{trackFiles[index]?.name || `Track ${index + 1}`}</Typography> */}
            </Box>
          ))}
        </Box>

        {/* Box 2: Selected tracks with move buttons */}
        <Typography mb={1}>
          Selected Tracks (Reorderable)
        </Typography>
        <Box
          sx={{
            display: "flex",
            flexDirection: "column",
            gap: 1,
            p: 2,
            border: "1px solid #ccc",
            borderRadius: "8px",
            backgroundColor: "#eef3ff",
          }}
        >
          {selectedTracks.map((track, index) => (
            <Box
              key={index}
              sx={{
                display: "flex",
                alignItems: "center",
                justifyContent: "space-between",
                gap: 2,
                padding: 1,
                backgroundColor: "#fff",
                borderRadius: "6px",
                boxShadow: "sm",
              }}
            >
              <Typography>{track.name}</Typography>
              {/* <Typography>{trackFiles[tracks.indexOf(track)]?.name || index}</Typography> */}
              <Box sx={{ display: "flex", gap: 1 }}>
                <IconButton onClick={() => handleMove(index, -1)}>
                  <ArrowUpwardIcon />
                </IconButton>
                <IconButton onClick={() => handleMove(index, 1)}>
                  <ArrowDownwardIcon />
                </IconButton>
              </Box>
            </Box>
          ))}
        </Box>

        <Box mt={4} display="flex" justifyContent="flex-end" gap={2}>
          <Button variant="outlined" onClick={onClose}>
            Cancel
          </Button>
          <Button variant="solid" onClick={handleConfirm}>
            Confirm Selection
          </Button>
        </Box>
      </Sheet>
    </Box>
  );
};

export default TrackSelector;
