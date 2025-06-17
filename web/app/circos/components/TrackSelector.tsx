import React, { useEffect, useState } from "react";
import {
  Box,
  Checkbox,
  Typography,
  IconButton,
  Sheet,
  Button,
  Divider,
} from "@mui/joy";
import { Collapse } from "@mui/material";
import ArrowDownwardIcon from "@mui/icons-material/ArrowDownward";
import ArrowForwardIcon from '@mui/icons-material/ArrowForward';
import FolderIcon from "@mui/icons-material/Folder";
import DescriptionIcon from "@mui/icons-material/Description";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import ExpandLessIcon from "@mui/icons-material/ExpandLess";
import { Track, TrackType } from "../config/track";
import { defaultAssemblyConfig, defaultChordConfig, defaultGlobalConfig, defaultLineConfig } from "../config/defaultConfigs";
import { useAppSelector } from "@/store/hooks";

interface FileEntry {
  name: string;
  type: "file" | "directory";
  size?: number;
  children?: FileEntry[];
}

interface TrackSelectorProps {
  onClose: () => void;
  onConfirm: (selectedTracks: Track[]) => void;
}

const isFasta = (track: Track) =>
  track.trackType === TrackType.Karyotype;

const isReferenceFasta = (track: Track) =>
  isFasta(track) && track.reference;

const TrackSelector: React.FC<TrackSelectorProps> = ({
  onClose,
  onConfirm,
}) => {
  const [files, setFiles] = useState<FileEntry[]>([]);
  const [selectedTracks, setSelectedTracks] = useState<Track[]>([]);
  const [expandedFolders, setExpandedFolders] = useState<Set<string>>(new Set());
  const space = useAppSelector((state) => state.space);
  const [trackMap, setTrackMap] = useState<Map<string, Track>>(new Map());

  useEffect(() => {

    const fetchFiles = async () => {
      try {
        const host = process.env.NEXT_PUBLIC_DJANGO_HOST;
        const res = await fetch(`${host}/api/timge/get_files_hierarchical/`, {
          method: "POST",
          headers: {
            "Content-Type": "application/json",
          },
          body: JSON.stringify({
            uuid: space.uuid,
            path: false,
          }),
        });
        const data = await res.json();
    
        if (data.status === "success") {
          const tracks = new Map<string, Track>();
    
          const walk = (entries: FileEntry[], parentPath = ""): FileEntry[] => {
            return entries.map((entry) => {
              const fullPath = parentPath ? `${parentPath}/${entry.name}` : entry.name;
    
              if (entry.type === "directory") {
                return {
                  ...entry,
                  children: entry.children ? walk(entry.children, fullPath) : [],
                };
              }
    
              const track = createTrackFromEntry(entry, fullPath);
              if (track) {
                tracks.set(entry.name, track);
              }
    
              return entry;
            });
          };
    
          const transformed = walk(data.entries);
          setFiles(transformed);
          setTrackMap(tracks);
        } else {
          console.error("Fetch error:", data.message);
        }
      } catch (err) {
        console.error("Request failed:", err);
      }
    };

    fetchFiles();
  }, [space]);

  // const toggleTrack = (track: Track) => {
  //   setSelectedTracks((prev) =>
  //     prev.includes(track)
  //       ? prev.filter((t) => t !== track)
  //       : [...prev, track]
  //   );
  // };

  const toggleTrack = (track: Track) => {
    setSelectedTracks((prev) => {
      const alreadySelected = prev.find((t) => t.name === track.name);

      if (alreadySelected) {
        return prev.filter((t) => t.name !== track.name);
      }

      const newSelection = [...prev, track];

      // Count how many FASTA tracks are selected
      const fastaTracks = newSelection.filter(isFasta);

      // If only one FASTA, mark it as reference
      if (fastaTracks.length === 1) {
        fastaTracks[0].reference = true;
      } else if (fastaTracks.length > 1) {
        // If more than one, ensure only one is reference
        // Keep the first one as reference, others not
        fastaTracks.forEach((t, i) => {
          t.reference = i === 0;
        });
      }

    return newSelection;
    });
  };

  const createTrackFromEntry = (entry: FileEntry, path: string): Track | null => {
    const lower = entry.name.toLowerCase();
  
    if (lower.endsWith(".bed") || lower.endsWith(".bedgraph")) {
      return {
        name: entry.name,
        trackType: TrackType.Line,
        config: { ...defaultLineConfig, filePath: path },
        data: {},
        path: path,
      };
    } else if (lower.endsWith(".fa") || lower.endsWith(".fasta")) {
      return {
        name: entry.name,
        trackType: TrackType.Karyotype,
        config: { ...defaultAssemblyConfig, filePath: path },
        data: {},
        path: path,
      };
    } else if (lower.endsWith(".bedpe")) {
      return {
        name: entry.name,
        trackType: TrackType.Chord,
        config: { ...defaultChordConfig, filePath: path },
        data: {},
        path: path,
      };
    }
  
    return null;
  };

  const toggleFolder = (path: string) => {
    const updated = new Set(expandedFolders);
    updated.has(path) ? updated.delete(path) : updated.add(path);
    setExpandedFolders(updated);
  };

  const renderTree = (entries: FileEntry[], parentPath = "") => {
    console.log("Rendering tree with entries:", entries);
    return entries.map((entry) => {
      const fullPath = parentPath ? `${parentPath}/${entry.name}` : entry.name;

      if (entry.type === "directory") {
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
              {entry.children && renderTree(entry.children, fullPath)}
            </Collapse>
          </Box>
        );
      }

      return (
        <Box key={fullPath} sx={{ display: "flex", alignItems: "center", ml: 4 }}>
        <Checkbox
          checked={selectedTracks.some((track) => track.name === entry.name)}
          disabled={!trackMap.has(entry.name)}
          onChange={() => {
            const track = trackMap.get(entry.name);
            if (track) toggleTrack(track);
          }}
        />
          <DescriptionIcon sx={{ml: 1}}/>
          <Typography sx={{ ml: 1 }}>{entry.name}</Typography>
        </Box>
      );
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

        <Typography mb={1}>Selected Tracks (Reorderable)</Typography>
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
              <Box sx={{ display: "flex", gap: 1 }}>
                <IconButton onClick={() => handleMove(index, -1)}>
                  <ArrowForwardIcon />
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
