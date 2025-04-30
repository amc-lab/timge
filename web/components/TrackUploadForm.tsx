"use client";
import {
  Box,
  Button,
  Typography,
  Stack,
  Divider,
  Card,
  IconButton,
  Sheet,
} from "@mui/joy";
import { Collapse } from "@mui/material";
import React, { useState } from "react";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import ExpandLessIcon from "@mui/icons-material/ExpandLess";
import FolderIcon from "@mui/icons-material/Folder";
import DescriptionIcon from "@mui/icons-material/Description";
import DataTrackFileUploadBox from "./TrackFileUpload";
import { useAppSelector } from "@/store/hooks";
import { useAppDispatch } from "@/store/hooks";
import { setLoading } from "@/store/features/site/siteSlice";
import DeleteIcon from '@mui/icons-material/Delete';

function TrackUploadForm({
  isOpen,
  onClose,
  tracks = [],
  onTrackUpload,
  onDeleteTrack,
}) {
  const _files = useAppSelector((state) => state.files);
  const dispatch = useAppDispatch();
  const [expandedFolders, setExpandedFolders] = useState(new Set());

  if (!isOpen) return null;

  const toggleFolder = (path) => {
    const updated = new Set(expandedFolders);
    updated.has(path) ? updated.delete(path) : updated.add(path);
    setExpandedFolders(updated);
  };

  const getFullPath = (track, parentPath = "") => {
    return parentPath ? `${parentPath}/${track.name}` : track.name;
  };

  const renderTree = (entries, parentPath = "") => {
    return entries.map((entry) => {
      const fullPath = getFullPath(entry, parentPath);
      if (entry.children) {
        return (
          <Box key={fullPath} sx={{ ml: 2 }}>
            <Box
              sx={{ display: "flex", alignItems: "center", cursor: "pointer" }}
              onClick={() => toggleFolder(fullPath)}
            >
              <IconButton size="sm">
                {expandedFolders.has(fullPath) ? <ExpandLessIcon /> : <ExpandMoreIcon />}
              </IconButton>
              <FolderIcon />
              <Typography sx={{ ml: 1 }}>{entry.name}</Typography>
            </Box>
            <Collapse in={expandedFolders.has(fullPath)}>
              {renderTree(entry.children, fullPath)}
            </Collapse>
          </Box>
        );
      }
      return (
        <Box
          key={fullPath}
          sx={{ display: "flex", justifyContent: "space-between", alignItems: "center", backgroundColor: "#f5f5f5", borderRadius: "8px", padding: "8px 12px", ml: 4, mt: 1 }}
        >
          <Box sx={{ display: "flex", alignItems: "center" }}>
            <DescriptionIcon />
            <Typography sx={{ ml: 1 }}>{entry.name}</Typography>
          </Box>
            <IconButton size="sm" color="danger" onClick={() => onDeleteTrack(entry.path)}>
              <DeleteIcon />
            </IconButton>
        </Box>
      );
    });
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
          sx={{ padding: 2, backgroundColor: "#f5f5f5", borderRadius: "8px", mb: 2 }}
        >
          <Typography level="body-md">
            Upload your own tracks here.
          </Typography>
        </Card>

        {_files.length > 0 && (
          <Box sx={{ mt: 3, mb: 3 }}>
            <Typography level="title-md" sx={{ mb: 1 }}>
              Uploaded Tracks
            </Typography>
            <Sheet sx={{ p: 2, border: "1px solid #ccc", borderRadius: "8px", backgroundColor: "#f9f9f9" }}>
              {renderTree(_files)}
            </Sheet>
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
          sx={{ padding: 2, backgroundColor: "#f5f5f5", borderRadius: "8px", mb: 2 }}
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
