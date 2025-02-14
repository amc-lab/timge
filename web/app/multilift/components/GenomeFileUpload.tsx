import React from "react";
import { useDropzone } from "react-dropzone";
import { Card, Typography } from "@mui/material";
import CloudUploadIcon from "@mui/icons-material/CloudUpload";

interface GenomeFileUploadBoxProps {
    genome: string;
    onGenomeFileUpload: (file: File) => void;
    onGenomeFileRemove: () => void;
}

const GenomeFileUploadBox: React.FC<GenomeFileUploadBoxProps> = ({ genome, onGenomeFileUpload, onGenomeFileRemove }) => {

    const { getRootProps, getInputProps, isDragActive } = useDropzone({
        onDrop: (acceptedFiles) => {
            if (acceptedFiles.length > 0) {
                setIsUploaded(true);
                setFileName(acceptedFiles[0].name);
                onGenomeFileUpload(acceptedFiles[0]);
            }
        },
    });

    const removeFile = () => {
        setIsUploaded(false);
        setFileName("");
        onGenomeFileRemove();
    }

    const [isUploaded, setIsUploaded] = React.useState(false);
    const [fileName, setFileName] = React.useState("");

    return (
        <Card
        sx={{
            width: "calc(25% - 1em)",
            aspectRatio: "1 / 1",
            display: "flex",
            alignItems: "center",
            justifyContent: "center",
            outline: "3px dashed #ccc",
            boxShadow: "none",
            flexDirection: "column",
            backgroundColor: isDragActive ? "#f0f0f0" : "#f9f9f9",
            // backgroundColor: "#f9f9f9",
            padding: "1em",
            textAlign: "center",
            "&:hover": {
                backgroundColor: "#f0f0f0",
            },
            marginTop: "0.5em",
            marginBottom: "0.5em",
            ...(!isUploaded && {
                cursor: "pointer",
            }),
            transition: "background-color 0.2s ease-in-out", // Smooth transition
        }}
        {...getRootProps()}
        >
            {
                !isUploaded ? (
                    <>
                        <input {...getInputProps()} type="file" />
                        <CloudUploadIcon sx={{ fontSize: 40, marginBottom: "0.5em", color: "#777" }} />
                        <Typography
                            sx={{
                            // fontSize: "0.8em",
                            }}
                        >
                            {`Drag & drop file for genome: ${genome}, or click to select`}
                        </Typography>
                    </>
                ) : (
                    <>
                        <Typography
                            sx={{
                            // fontSize: "0.8em",
                            }}
                        >
                            {`Uploaded genome: ${fileName}`}
                        </Typography>
                        <Typography
                            sx={{
                                marginTop: "0.5em",
                                textDecoration: "underline",
                                color: "blue",
                                cursor: "pointer",
                            }}
                            onClick={removeFile}
                        >
                            Remove file
                        </Typography>
                    </>
                )
            }

        </Card>
    );
};

export default GenomeFileUploadBox;
