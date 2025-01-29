import React, { useState, ChangeEvent } from "react";
import { Card } from "@mui/joy";
import Button from "@mui/joy/Button";
import Select from "@mui/joy/Select";
import Option from "@mui/joy/Option";

interface FileEntry {
  name: string;
  data: any;
  trackType: string;
}

interface FileUploadFormProps {
  onSubmit: (files: FileEntry[]) => void;
}

const FileUploadForm: React.FC<FileUploadFormProps> = ({ onSubmit }) => {
  const [files, setFiles] = useState<FileEntry[]>([]);

  const trackTypes = ["Karyotype", "Bar", "Chord"];

  const handleFileUpload = (file: File) => {
    const reader = new FileReader();
    reader.onload = (event) => {
      try {
        const json = JSON.parse(event.target?.result as string);
        const newFile: FileEntry = {
          name: file.name,
          data: json,
          trackType: trackTypes[0],
        };
        setFiles((prevFiles) => [...prevFiles, newFile]);
      } catch (error) {
        console.error("Invalid JSON file:", error);
        alert("Uploaded file is not valid JSON.");
      }
    };
    reader.readAsText(file);
  };

  const handleFileChange = (event: ChangeEvent<HTMLInputElement>) => {
    if (event.target.files) {
      Array.from(event.target.files).forEach((file) => handleFileUpload(file));
    }
  };

  const handleTrackTypeChange = (index: number, newTrackType: string) => {
    setFiles((prevFiles) => {
      const updatedFiles = [...prevFiles];
      updatedFiles[index].trackType = newTrackType;
      return updatedFiles;
    });
  };

  const handleFileRemove = (index: number) => {
    setFiles((prevFiles) => prevFiles.filter((_, i) => i !== index));
  };

  const moveFile = (fromIndex: number, toIndex: number) => {
    setFiles((prevFiles) => {
      const updatedFiles = [...prevFiles];
      const [movedFile] = updatedFiles.splice(fromIndex, 1);
      updatedFiles.splice(toIndex, 0, movedFile);
      return updatedFiles;
    });
  };

  const handleSubmit = () => {
    onSubmit(files);
  };

  return (
    <Card>
      <h2>Upload and Configure Files</h2>

      <input type="file" multiple onChange={handleFileChange} accept=".json" />

      <div style={{ marginTop: "20px" }}>
        {files.map((fileEntry, index) => (
          <div
            key={index}
            style={{
              display: "flex",
              alignItems: "center",
              marginBottom: "10px",
            }}
          >
            <span style={{ marginRight: "10px" }}>{fileEntry.name}</span>

            <Select
              placeholder="Select visualisation type"
              value={fileEntry.trackType}
              onChange={(event, newValue) =>
                newValue && handleTrackTypeChange(index, newValue)
              }
              style={{ marginRight: "10px" }}
            >
              {trackTypes.map((type, idx) => (
                <Option key={idx} value={type}>
                  {type}
                </Option>
              ))}
            </Select>

            <Button
              onClick={() => handleFileRemove(index)}
              style={{ backgroundColor: "red", color: "white", border: "none" }}
            >
              Remove
            </Button>

            {index > 0 && (
              <Button
                onClick={() => moveFile(index, index - 1)}
                style={{ marginLeft: "10px" }}
              >
                Move Up
              </Button>
            )}

            {index < files.length - 1 && (
              <Button
                onClick={() => moveFile(index, index + 1)}
                style={{ marginLeft: "10px" }}
              >
                Move Down
              </Button>
            )}
          </div>
        ))}
      </div>

      {files.length > 0 && <Button onClick={handleSubmit}>Submit</Button>}
    </Card>
  );
};

export default FileUploadForm;
