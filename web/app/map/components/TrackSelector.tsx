import React, { useEffect, useState } from "react";
import {
  Box,
  Typography,
  IconButton,
} from "@mui/joy";
import { Collapse } from "@mui/material";
import FolderIcon from "@mui/icons-material/Folder";
import DescriptionIcon from "@mui/icons-material/Description";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import ExpandLessIcon from "@mui/icons-material/ExpandLess";
import { useAppSelector } from "@/store/hooks";
import TrackSelectorModal from "@/components/TrackSelectorModal";
import { API_BASE_URL } from "@/app/config/env";

interface FileEntry {
  name: string;
  type: "file" | "directory";
  size?: number;
  children?: FileEntry[];
}

interface TrackSelectorProps {
  onClose: () => void;
  onConfirm: (referencePath: string | null, trackPath: string | null) => void;
}

const TrackSelector: React.FC<TrackSelectorProps> = ({ onClose, onConfirm }) => {
  const [files, setFiles] = useState<FileEntry[]>([]);
  const [expandedFolders, setExpandedFolders] = useState<Set<string>>(new Set());
  const [referencePath, setReferencePath] = useState<string | null>(null);
  const [trackPath, setTrackPath] = useState<string | null>(null);
  const space = useAppSelector((state) => state.space);

  useEffect(() => {
    const fetchFiles = async () => {
      try {
        const res = await fetch(`${API_BASE_URL}/api/timge/get_files_hierarchical/`, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ uuid: space.uuid, path: false }),
        });
        const data = await res.json();

        if (data.status === "success") {
          const walk = (entries: FileEntry[], parentPath = ""): FileEntry[] =>
            entries.map((entry) => {
              const fullPath = parentPath ? `${parentPath}/${entry.name}` : entry.name;
              if (entry.type === "directory") {
                return {
                  ...entry,
                  children: entry.children ? walk(entry.children, fullPath) : [],
                };
              }
              return entry;
            });

          setFiles(walk(data.entries));
        } else {
          console.error("Fetch error:", data.message);
        }
      } catch (err) {
        console.error("Request failed:", err);
      }
    };

    fetchFiles();
  }, [space]);

  const toggleFolder = (path: string) => {
    const updated = new Set(expandedFolders);
    updated.has(path) ? updated.delete(path) : updated.add(path);
    setExpandedFolders(updated);
  };

  const handleFileSelect = (path: string, name: string) => {
    const lower = name.toLowerCase();
    if (lower.endsWith(".fa") || lower.endsWith(".fasta")) {
      setReferencePath(path);
    } else {
      setTrackPath(path);
    }
  };

  const renderTree = (entries: FileEntry[], parentPath = "") =>
    entries.map((entry) => {
      const fullPath = parentPath ? `${parentPath}/${entry.name}` : entry.name;

      if (entry.type === "directory") {
        return (
          <Box key={fullPath} sx={{ ml: 2 }}>
            <Box sx={{ display: "flex", alignItems: "center", cursor: "pointer" }} onClick={() => toggleFolder(fullPath)}>
              <IconButton size="sm">
                {expandedFolders.has(fullPath) ? <ExpandLessIcon /> : <ExpandMoreIcon />}
              </IconButton>
              <FolderIcon />
              <Typography sx={{ ml: 1 }}>{entry.name}</Typography>
            </Box>
            <Collapse in={expandedFolders.has(fullPath)}>
              {entry.children && renderTree(entry.children, fullPath)}
            </Collapse>
          </Box>
        );
      }

      return (
        <Box
          key={fullPath}
          sx={{
            display: "flex",
            alignItems: "center",
            ml: 4,
            cursor: "pointer",
            backgroundColor: fullPath === referencePath || fullPath === trackPath ? "#e0f7fa" : "transparent",
            p: 1,
            borderRadius: 1,
          }}
          onClick={() => handleFileSelect(fullPath, entry.name)}
        >
          <DescriptionIcon />
          <Typography sx={{ ml: 1 }}>{entry.name}</Typography>
        </Box>
      );
    });

  return (
    <TrackSelectorModal
      title="Select Reference and Track File"
      description="Choose a reference FASTA and corresponding track file."
      onClose={onClose}
      onConfirm={() => onConfirm(referencePath, trackPath)}
      confirmDisabled={!referencePath || !trackPath}
    >
      <Box
        sx={{
          border: "1px solid #ccc",
          borderRadius: "8px",
          backgroundColor: "#f9f9f9",
          p: 2,
          mb: 3,
        }}
      >
        {renderTree(files)}
      </Box>

      <Typography>
        Reference File: {referencePath || "None selected"}
      </Typography>
      <Typography>Track File: {trackPath || "None selected"}</Typography>
    </TrackSelectorModal>
  );
};

export default TrackSelector;
