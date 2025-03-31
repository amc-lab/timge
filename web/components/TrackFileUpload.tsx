import React from "react";
import { useDropzone } from "react-dropzone";
import { Card, Typography } from "@mui/material";
import CloudUploadIcon from "@mui/icons-material/CloudUpload";

interface DataTrackFileUploadBoxProps {
    onDataTrackFileUpload: (fileList: FileList) => void;
}

const DataTrackFileUploadBox: React.FC<DataTrackFileUploadBoxProps> = ({ onDataTrackFileUpload }) => {

    const { getRootProps, getInputProps, isDragActive } = useDropzone({
        onDrop: (acceptedFiles) => {
            if (acceptedFiles.length > 0) {
                const fileList = new DataTransfer();
                acceptedFiles.forEach(file => fileList.items.add(file));
                onDataTrackFileUpload(fileList.files);
            }
        },
    });

    return (
        <Card
        sx={{
            width: "100%",
            display: "flex",
            aspectRatio: "5 / 1",
            alignItems: "center",
            justifyContent: "center",
            outline: "2px dashed #ccc",
            boxShadow: "none",
            flexDirection: "column",
            backgroundColor: isDragActive ? "#f0f0f0" : "#f9f9f9",
            padding: "2em",
            textAlign: "center",
            "&:hover": {
                backgroundColor: "#f0f0f0",
            },
            marginTop: "0.5em",
            marginBottom: "0.5em",
            cursor: "pointer",
            transition: "background-color 0.2s ease-in-out", // Smooth transition
        }}
        {...getRootProps()}
        >
            <>
                <input {...getInputProps()} type="file" />
                <CloudUploadIcon sx={{ fontSize: 40, marginBottom: "0.5em", color: "#777" }} />
                <Typography
                >
                    {`Drag & drop track file(s)`}
                </Typography>
            </>
        </Card>
    );
};

export default DataTrackFileUploadBox;
