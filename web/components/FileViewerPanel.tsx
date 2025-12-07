"use client";
import React, { useEffect, useState } from 'react';
import { useDropzone } from "react-dropzone";
import { Card, LinearProgress } from "@mui/material";
import { Box, Typography, IconButton } from '@mui/joy';
import CloseIcon from '@mui/icons-material/Close';
import FolderIcon from '@mui/icons-material/Folder';
import InsertDriveFileIcon from '@mui/icons-material/InsertDriveFile';

import { useAppDispatch, useAppSelector } from "@/store/hooks";
import { setSidebarExpanded, setWorkingDirectory } from '@/store/features/space/spaceSlice';
import { fetchFiles } from '@/store/features/files/fileSlice';
import { uploadTrackFiles } from '@/app/utils/fileUtils';
import { API_BASE_URL } from "@/app/config/env";

import ContextMenu, { ContextMenuOption } from './FileViewerPanelContextMenu';

interface FileEntry {
  name: string;
  isDirectory: boolean;
  size?: number;
  children?: FileEntry[];
  path: string;
}

const FileViewerPanel: React.FC<{
  refreshView: boolean;
  setRefreshView: (b: boolean) => void;
}> = ({ refreshView, setRefreshView }) => {
  const dispatch = useAppDispatch();
  const space = useAppSelector((s) => s.space);
  const _files = useAppSelector((s) => s.files);

  const [workingDir, setWorkingDir] = useState<FileEntry[]>([]);
  const [isUploading, setIsUploading] = useState(true);

  const [menuState, setMenuState] = useState<{
    visible: boolean;
    x: number;
    y: number;
    type: "background" | "file";
    file?: FileEntry;
  }>({
    visible: false,
    x: 0,
    y: 0,
    type: "background",
  });

  // 1) build the correct subdirectory
  const getDirectoryInWorkingDir = (): FileEntry[] => {
    let dir = _files as FileEntry[];
    for (const part of space.config.working_directory) {
      const found = dir.find((f) => f.name === part && f.children);
      if (found && found.children) dir = found.children;
      else return _files as FileEntry[];
    }
    return dir;
  };

  useEffect(() => {
    setWorkingDir(getDirectoryInWorkingDir());
    setIsUploading(false);
  }, [_files, space.config.working_directory]);

  const { getRootProps, getInputProps, isDragActive } = useDropzone({
    onDrop: async (files) => {
      if (!files.length) return;
      setIsUploading(true);
      try {
        await uploadTrackFiles(space, files);
        dispatch(fetchFiles({ uuid: space.uuid, path: [] }));
      } catch (err) {
        console.error(err);
      } finally {
        setIsUploading(false);
      }
    },
    noClick: true,
    noKeyboard: true,
  });

  const handleContextMenu = (
    e: React.MouseEvent,
    type: "background" | "file",
    file?: FileEntry
  ) => {
    e.preventDefault();
    e.stopPropagation();
    setMenuState({
      visible: true,
      x: e.clientX,
      y: e.clientY,
      type,
      file,
    });
  };
  const handleCloseMenu = () =>
    setMenuState((ms) => ({ ...ms, visible: false }));

const handleDownload = async (file: FileEntry) => {
    const url = new URL(`${API_BASE_URL}/api/timge/download/`);
    url.searchParams.append("uuid", space.uuid);
    url.searchParams.append("path", file.path);

    const a = document.createElement("a");
    a.href = url.toString();
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
};

  const backgroundOptions: ContextMenuOption[] = [
    {
      label: "Refresh",
      onClick: () => dispatch(fetchFiles({ uuid: space.uuid, path: space.config.working_directory })),
    },
  ];

  const fileOptions: ContextMenuOption[] = [
    {
        label: "Download",
        onClick: () => {
          if (menuState.file) {
            handleDownload(menuState.file).catch((e) => {
              console.error(e);
              alert(`Download error: ${e.message}`);
            });
          }
        },
      },
  ];

  return (
    <Card
      component="div"
      {...getRootProps()}
      onContextMenu={(e) => handleContextMenu(e, "background")}
      sx={{
        width: "17.5%",
        minHeight: "100vh",
        borderRadius: 0,
        display: "flex",
        flexDirection: "column",
        position: "relative",
        backgroundColor: isDragActive ? "#f4f4f4" : "#FBFCFE",
        outline: isDragActive ? "2px solid #1e90ff" : "none",
        transition: "background-color 0.2s ease-in-out",
        p: 2,
      }}
    >
      <ContextMenu
        visible={menuState.visible}
        x={menuState.x}
        y={menuState.y}
        onClose={handleCloseMenu}
        options={
          menuState.type === "background"
            ? backgroundOptions
            : fileOptions
        }
      />

      {isUploading && (
        <LinearProgress
          sx={{
            position: "absolute",
            top: 0,
            left: 0,
            right: 0,
            zIndex: 1,
          }}
        />
      )}

      <input {...getInputProps()} type="file" />

      <Box sx={{ display: "flex", alignItems: "center", justifyContent: "space-between", pb: 1, borderBottom: 1, borderColor: "divider" }}>
        <Typography>File Viewer</Typography>
        <IconButton size="sm" variant="plain" onClick={() => dispatch(setSidebarExpanded(false))}>
          <CloseIcon />
        </IconButton>
      </Box>

      {space.config.working_directory.length > 0 && (
        <Box
          onClick={() =>
            dispatch(setWorkingDirectory(space.config.working_directory.slice(0, -1)))
          }
          sx={{
            display: "flex",
            alignItems: "center",
            py: 1,
            borderBottom: 1,
            borderColor: "divider",
            cursor: "pointer",
            "&:hover": { bgcolor: "background.level1" },
            justifyContent: "center",
          }}
        >
          <Typography>…</Typography>
        </Box>
      )}

      <Box sx={{ overflowY: "auto", flexGrow: 1 }}>
        {workingDir.map((file, i) => (
          <Box
            key={i}
            onClick={() => {
              if (file.isDirectory) {
                dispatch(
                  setWorkingDirectory([...space.config.working_directory, file.name])
                );
              }
            }}
            onContextMenu={(e) => handleContextMenu(e, "file", file)}
            sx={{
                display: "flex",
                alignItems: "center",
                paddingY: "1em",
                borderBottom: "1px solid #ccc",
                fontSize: "0.8em",
                ":hover": {
                    backgroundColor: "#ECF0F1",
                },
                cursor: "pointer",
            }}
          >
            {file.isDirectory ? (
              <FolderIcon fontSize="small" sx={{ mr: 1 }} />
            ) : (
              <InsertDriveFileIcon fontSize="small" sx={{ mr: 1 }} />
            )}
            <Box sx={{
                width: "60%",
            }}>
            <Typography sx={{ 
                flexGrow: 1, 
                textOverflow: "ellipsis", 
                overflow: "hidden", 
                fontSize: "1em",
                whiteSpace: "nowrap",
            }}>
              {file.name}
            </Typography>
            </Box>
            <Box
                sx={{
                    marginLeft: "auto",
                    display: "flex",
                    alignItems: "center",
                }}
            >
            {!file.isDirectory && file.size != null && (
              <Typography sx={{ ml: 1, fontSize: "0.9em" }}>
                {(file.size / 1e6).toFixed(2)} MB
              </Typography>
            )}
            </Box>
          </Box>
        ))}
      </Box>
    </Card>
  );
};

export default FileViewerPanel;
