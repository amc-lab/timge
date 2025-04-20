"use client";
import React, { useEffect, useState } from 'react';
import { useDropzone } from "react-dropzone";
import {Card, LinearProgress} from "@mui/material";
import { Box, Typography, IconButton } from '@mui/joy';
import CloseIcon from '@mui/icons-material/Close';
import { useAppDispatch } from "@/store/hooks";
import { useAppSelector } from "@/store/hooks";
import { setSidebarExpanded } from '@/store/features/space/spaceSlice';
import { uploadTrackFiles } from '@/app/utils/fileUtils';

interface FileViewerPanelProps {
    files: any[];
    triggerFileRefresh: () => void;
}

const FileViewerPanel = ({ files, triggerFileRefresh }: FileViewerPanelProps) => {
    const dispatch = useAppDispatch();
    const space = useAppSelector((state) => state.space);

    const handleFileDrop = async (acceptedFiles: File[]) => {
        setIsUploading(true);
        if (acceptedFiles.length > 0) {
        await uploadTrackFiles(space, Array.from(acceptedFiles));
        await new Promise(resolve => setTimeout(resolve, 1000));
        triggerFileRefresh();
        }
    };

    const { getRootProps, getInputProps, isDragActive } = useDropzone({
        onDrop: handleFileDrop,
        noClick: true,
        noKeyboard: true,
    });

    const [isUploading, setIsUploading] = useState(true);

    useEffect(() => {
        setIsUploading(false);
        console.log(files);
    }, [files]);

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
    >
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

                {files.map((file, index) => (
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
                    >
                        <Box sx={{ width: "60%" }}>
                            <Typography
                                sx={{
                                    // color: "black",
                                    fontSize: "1em",
                                    marginLeft: "5px",
                                    textOverflow: "ellipsis",
                                    whiteSpace: "nowrap",
                                    overflow: "hidden",
                                }}
                            >
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
                            <Typography
                                sx={{
                                    color: "black",
                                    fontSize: "1em",
                                    marginLeft: "10px",
                                }}
                            >
                                {(file.size / (1024 * 1024)).toFixed(2)} MB
                            </Typography>
                        </Box>
                    </Box>
                ))}
            </Box>
        </Card>
    );
};

export default FileViewerPanel;
