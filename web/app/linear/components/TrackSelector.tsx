import React, { useEffect, useState } from "react";
import {
  Box,
  Typography,
  IconButton,
  Sheet,
  Button,
  Divider,
  Checkbox,
} from "@mui/joy";
import { Collapse } from "@mui/material";
import FolderIcon from "@mui/icons-material/Folder";
import DescriptionIcon from "@mui/icons-material/Description";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import ExpandLessIcon from "@mui/icons-material/ExpandLess";
import { useAppSelector } from "@/store/hooks";

interface FileEntry {
  name: string;
  type: "file" | "directory";
  size?: number;
  children?: FileEntry[];
}

interface TrackSelectorProps {
  onClose: () => void;
  onConfirm: (referencePath: string, trackPaths: string[]) => void;
}

const TrackSelector: React.FC<TrackSelectorProps> = ({ onClose, onConfirm }) => {
  const [files, setFiles] = useState<FileEntry[]>([]);
  const [expandedFolders, setExpandedFolders] = useState<Set<string>>(new Set());
  const [referencePath, setReferencePath] = useState<string | null>(null);
  const [selectedTracks, setSelectedTracks] = useState<Set<string>>(new Set());
  const space = useAppSelector((state) => state.space);

  useEffect(() => {
    const fetchFiles = async () => {
      try {
        const host = process.env.NEXT_PUBLIC_DJANGO_HOST;
        const res = await fetch(`${host}/api/timge/get_files_hierarchical/`, {
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

  const handleReferenceSelect = (path: string) => {
    setReferencePath(path);
  };

  const handleTrackToggle = (path: string) => {
    const updated = new Set(selectedTracks);
    updated.has(path) ? updated.delete(path) : updated.add(path);
    setSelectedTracks(updated);
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

      const isReference = entry.name.endsWith(".fa") || entry.name.endsWith(".fasta");

      return (
        <Box
          key={fullPath}
          sx={{
            display: "flex",
            alignItems: "center",
            ml: 4,
            cursor: isReference ? "pointer" : "default",
            backgroundColor:
              (fullPath === referencePath && isReference) || selectedTracks.has(fullPath)
                ? "#e0f7fa"
                : "transparent",
            p: 1,
            borderRadius: 1,
          }}
          onClick={() => isReference && handleReferenceSelect(fullPath)}
        >
          {isReference ? (
            <>
              <DescriptionIcon />
              <Typography sx={{ ml: 1 }}>{entry.name} (Reference)</Typography>
            </>
          ) : (
            <>
              <Checkbox
                size="sm"
                checked={selectedTracks.has(fullPath)}
                onChange={() => handleTrackToggle(fullPath)}
              />
              <DescriptionIcon sx={{ml : 1}} />
              <Typography sx={{ ml: 0.5 }}>{entry.name}</Typography>
            </>
          )}
        </Box>
      );
    });

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
          Select Reference Genome and Tracks
        </Typography>
        <Divider sx={{ mb: 2 }} />

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

        <Typography>Reference File: {referencePath || "None selected"}</Typography>
        <Typography>Tracks Selected: {selectedTracks.size}</Typography>

        <Box mt={4} display="flex" justifyContent="flex-end" gap={2}>
          <Button variant="outlined" onClick={onClose}>
            Cancel
          </Button>
          <Button
            variant="solid"
            onClick={() => {
              if (referencePath) onConfirm(referencePath, Array.from(selectedTracks));
            }}
            disabled={!referencePath}
          >
            Confirm Selection
          </Button>
        </Box>
      </Sheet>
    </Box>
  );
};

export default TrackSelector;