import {
  Box,
  Button,
  Typography,
  Stack,
  Divider,
  Card,
} from "@mui/joy";
import React from "react";
import DataTrackFileUploadBox from "./TrackFileUpload";

function TrackUploadForm({
  isOpen,
  onClose,
  tracks = [],
  onTrackUpload,
  onDeleteTrack,
}) {
  if (!isOpen) return null;

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
        alignItems: "center",
        justifyContent: "center",
        zIndex: 1300,
      }}
    >
      <Box
        sx={{
          position: "relative",
          backgroundColor: "#fff",
          borderRadius: "12px",
          padding: 4,
          minWidth: "50vw",
          boxShadow: "lg",
          maxHeight: "80vh",
          overflowY: "auto",
        }}
      >
        <Button
          onClick={onClose}
          size="sm"
          variant="plain"
          color="neutral"
          sx={{
            position: "absolute",
            top: 12,
            right: 12,
            fontSize: "1.25rem",
            minWidth: "unset",
            padding: "4px 8px",
            lineHeight: 1,
          }}
        >
          &times;
        </Button>

        <Typography level="title-lg" sx={{ mb: 1 }}>
          My Tracks
        </Typography>
        <Divider sx={{ mb: 2 }} />

        <Card
          variant="outlined"
          sx={{
            padding: 2,
            backgroundColor: "#f5f5f5",
            borderRadius: "8px",
            mb: 2,
          }}
        >
          <Typography level="body-md">
            Upload your own tracks here.
          </Typography>
        </Card>

        {tracks.length > 0 && (
          <Box sx={{ mt: 3, mb: 3 }}>
            <Typography level="title-md" sx={{ mb: 1 }}>
              Uploaded Tracks
            </Typography>
            <Stack spacing={1}>
              {tracks.map((track, index) => (
                <Box
                  key={index}
                  sx={{
                    display: "flex",
                    justifyContent: "space-between",
                    alignItems: "center",
                    backgroundColor: "#f5f5f5",
                    borderRadius: "8px",
                    padding: "8px 12px",
                  }}
                >
                  <Typography>{track.name}</Typography>
                  <Button
                    size="sm"
                    color="danger"
                    onClick={() => onDeleteTrack(track.name)}
                  >
                    Delete
                  </Button>
                </Box>
              ))}
            </Stack>
          </Box>
        )}

        <DataTrackFileUploadBox
          onDataTrackFileUpload={(fileList) => {
            const files = Array.from(fileList);
            onTrackUpload(files);
          }}
        />

        <Typography level="title-lg" sx={{ mt: 2, mb: 1 }}>
          Preset Tracks
        </Typography>
        <Divider sx={{ mb: 2 }} />
        <Card
          variant="outlined"
          sx={{
            padding: 2,
            backgroundColor: "#f5f5f5",
            borderRadius: "8px",
            mb: 2,
          }}
        >
          <Typography level="body-md">
            Choose from the provided preset tracks.
          </Typography>
        </Card>
      </Box>
    </Box>
  );
}

export default TrackUploadForm;
