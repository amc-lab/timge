"use client";
import React, { useEffect, useState } from 'react';
import { useDropzone } from "react-dropzone";
import {Card, LinearProgress} from "@mui/material";
import { Box, Typography, IconButton, Breadcrumbs } from '@mui/joy';
import CloseIcon from '@mui/icons-material/Close';
import { useAppDispatch } from "@/store/hooks";
import { useAppSelector } from "@/store/hooks";
import { setSidebarExpanded, setWorkingDirectory } from '@/store/features/space/spaceSlice';
import { getTrackFiles, uploadTrackFiles } from '@/app/utils/fileUtils';
import FolderIcon from '@mui/icons-material/Folder';
import InsertDriveFileIcon from '@mui/icons-material/InsertDriveFile';
import ContextMenu from './FileViewerPanelContextMenu';
import { fetchFiles } from '@/store/features/files/fileSlice';

const FileViewerPanel = ({refreshView, setRefreshView}) => {
    const dispatch = useAppDispatch();
    const space = useAppSelector((state) => state.space);
    const _files = useAppSelector((state) => state.files);
    const [workingDir, setWorkingDir] = useState([]);
    const handleFileDrop = async (acceptedFiles: File[]) => {
        setIsUploading(true);
        if (acceptedFiles.length > 0) {
            await uploadTrackFiles(space, Array.from(acceptedFiles));
            await new Promise(resolve => setTimeout(resolve, 1000));
        }
        dispatch(fetchFiles({uuid: space.uuid, path: []}));
        setIsUploading(false);
    };

    const getDirectoryInWorkingDir = () => {
        let _dir = _files;
        for (let i = 0; i < space.config.working_directory.length; i++) {
            const foundDir = _dir.find((file) => file.name === space.config.working_directory[i]);
            if (foundDir && foundDir.children) {
                _dir = foundDir.children;
            } else {
                return _files;
            }
        }
        return _dir;
    }

    useEffect(() => {
        setWorkingDir(getDirectoryInWorkingDir());
        setIsUploading(false);
    }
    , [_files, space.config.working_directory]);

    const { getRootProps, getInputProps, isDragActive } = useDropzone({
        onDrop: handleFileDrop,
        noClick: true,
        noKeyboard: true,
    });

    const [isUploading, setIsUploading] = useState(true);
    const [menuState, setMenuState] = useState({
        visible: false,
        x: 0,
        y: 0,
      });

      const handleRightClick = (e: React.MouseEvent) => {
        e.preventDefault();
        setMenuState({
          visible: true,
          x: e.clientX,
          y: e.clientY,
        });
      };
    
      const handleClose = () => {
        setMenuState((prev) => ({ ...prev, visible: false }));
      };

    return (
    <Card
        component="div"
        {...getRootProps()}
        sx={{
            width: "17.5%",
            minHeight: "100vh",
            borderRadius: "0",
            display: "flex",
            flexDirection: "column",
            position: "relative",
            backgroundColor: isDragActive ? "#f4f4f4" : "#FBFCFE",
            outline: isDragActive ? "2px solid #1e90ff" : "none",
            transition: "background-color 0.2s ease-in-out",
            padding: "1em",
        }}
        onContextMenu={handleRightClick}
    >
        <ContextMenu
            visible={menuState.visible}
            x={menuState.x}
            y={menuState.y}
            onClose={handleClose}
            onCreateFolder={() => {
                console.log("Create folder clicked");
            }}
            onDelete={() => {
                console.log("Delete clicked");
            }}
            onRename={() => {
                console.log("Rename clicked");
            }}
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
            />)}
        <input {...getInputProps()} />
                
            <Box
                sx={{
                    overflowY: "auto",
                    flexGrow: 1,
                }}
            >
                <Box
                    sx={{
                        display: "flex",
                        alignItems: "center",
                        justifyContent: "space-between",
                        paddingY: "0.5em",
                        fontSize: "0.8em",
                        borderBottom: "1px solid #ccc",
                    }}
                >
                    <Typography
                        sx={{
                            color: "black",
                            fontSize: "1.2em",
                            fontWeight: "bold",
                            marginLeft: "5px",
                        }}
                    >
                        File Viewer
                    </Typography>
                    <IconButton
                        variant="plain"
                        size="sm"
                        onClick={() => dispatch(setSidebarExpanded(false))}
                    >
                        <CloseIcon />
                    </IconButton>
                </Box>

                {/* <Breadcrumbs
                    aria-label="breadcrumb"
                    sx={{
                        marginTop: "1em",
                        marginBottom: "1em",
                        fontSize: "0.8em",
                    }}
                >
                    {directory.map((item: string) => (
                        <Typography key={item} sx={{ color: "black" }}>
                            {item}
                        </Typography>
                    ))}
                </Breadcrumbs> */}
                {space.config.working_directory.length > 0 && (
                <Box
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
                    onClick={() => {
                        dispatch(setWorkingDirectory(space.config.working_directory.slice(0, -1)));
                    }}
                >
                    <Box sx={{ 
                        display: "flex", 
                        alignItems: "center",
                        width: "100%",
                        justifyContent: "center",
                    }}>
                        <Typography
                            sx={{
                                fontSize: "1em",
                                textOverflow: "ellipsis",
                                whiteSpace: "nowrap",
                                overflow: "hidden",
                            }}
                        >
                            ...
                        </Typography>
                    </Box>
                </Box>
                )}

                {workingDir.map((file, index) => {
                    const isDir = file.isDirectory;
                    return (
                        <Box
                            key={index}
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
                            onClick={() => {
                                if (isDir) {
                                    dispatch(setWorkingDirectory([...space.config.working_directory, file.name]));
                                }
                            }}
                        >
                            <Box sx={{ marginRight: "10px", display: "flex", alignItems: "center" }}>
                                {isDir ? <FolderIcon fontSize='small' /> : <InsertDriveFileIcon fontSize='small' />}
                            </Box>

                            <Box sx={{ width: "60%" }}>
                                <Typography
                                    sx={{
                                        fontSize: "1em",
                                        textOverflow: "ellipsis",
                                        whiteSpace: "nowrap",
                                        overflow: "hidden",
                                    }}
                                >
                                    {file.name}
                                </Typography>
                            </Box>

                            {/* File size (only for files) */}
                            {!isDir && (
                                <Box
                                    sx={{
                                        marginLeft: "auto",
                                        display: "flex",
                                        alignItems: "center",
                                    }}
                                >
                                    <Typography
                                        sx={{
                                            color: "black",
                                            fontSize: "1em",
                                            marginLeft: "10px",
                                        }}
                                    >
                                        {(file.size / (1000 * 1000)).toFixed(2)} MB
                                    </Typography>
                                </Box>
            )}
        </Box>
    );
})}
            </Box>
        </Card>
    );
};

export default FileViewerPanel;
