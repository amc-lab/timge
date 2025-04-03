import React, { useState } from "react";
import {
  Box,
  Checkbox,
  Typography,
  Select,
  Option,
  IconButton,
  Sheet,
  Button,
  Divider,
  Card,
} from "@mui/joy";
import DragIndicatorIcon from "@mui/icons-material/DragIndicator";
import { Track } from "../config/track";

interface TrackSelectorProps {
  tracks: Track[];
  trackFiles: any[];
  onClose: () => void;
  onConfirm: (selectedTracks: Track[]) => void;
}

const visualizationOptions = ["Linear View", "Circular", "Heatmap"];

const TrackSelector: React.FC<TrackSelectorProps> = ({ tracks, trackFiles, onClose, onConfirm }) => {

  const handleMove = (index: number, direction: number) => {
    // const newIndex = index + direction;
    // if (newIndex < 0 || newIndex >= trackOrder.length) return;
    // const newOrder = [...trackOrder];
    // const [movedItem] = newOrder.splice(index, 1);
    // newOrder.splice(newIndex, 0, movedItem);
    // setTrackOrder(newOrder);
  };

  const [selectedTracks, setSelectedTracks] = useState<Set<Track>>(new Set());

  const handleConfirm = () => {
    onConfirm(selectedTracks ? Array.from(selectedTracks) : []);
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

        <Box sx={{ display: "flex", flexDirection: "column", gap: 1 }}>
          {tracks.map((track, index) => {

            return (
              <Box
                key={index}
                sx={{
                  display: "flex",
                  alignItems: "center",
                  gap: 1,
                  padding: 1,
                  backgroundColor: "#f5f5f5",
                  borderRadius: "8px",
                  outline: "1px solid #ccc",
                }}
              >
                <IconButton onClick={() => handleMove(index, -1)}>
                  ↑
                </IconButton>
                <IconButton onClick={() => handleMove(index, 1)}>
                  ↓
                </IconButton>
                <Checkbox
                  checked={selectedTracks.has(track)}
                  onChange={(e) =>
                    setSelectedTracks((prev) => {
                      const newSet = new Set(prev);
                      if (e.target.checked) {
                        newSet.add(track);
                      } else {
                        newSet.delete(track);
                      }
                      return newSet;
                    }
                  )}
                />
                <Typography sx={{ minWidth: 120 }}>{trackFiles[index].name}</Typography>
                {/* <Select
                  value={visualizationTypes[track.id]}
                  onChange={(_, value) =>
                    setVisualizationTypes((prev) => ({
                      ...prev,
                      [track.id]: value || visualizationOptions[0],
                    }))
                  }
                  sx={{ minWidth: 150 }}
                >
                  {visualizationOptions.map((opt, idx) => (
                    <Option key={idx} value={opt}>
                      {opt}
                    </Option>
                  ))}
                </Select> */}
              </Box>
            );
          })}
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
