import { createSlice, PayloadAction, createAsyncThunk } from "@reduxjs/toolkit";
import { File } from "./types";
import { API_BASE_URL } from "@/app/config/env";

const queryFiles = async (uuid: string, path: string[]): Promise<File[]> => {
    const res = await fetch(`${API_BASE_URL}/api/timge/get_files_hierarchical/`, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({
          uuid: uuid,
          path: path,
        }),
    })
    .then(response => response.json())
    .then(data => {
        if (data.status === "success") {
            const formatEntries = (entries: any[], currentPath: string[]): File[] => {
                return entries.map((entry) => {
                    const entryPath = [...currentPath, entry.name];
                    const file: File = {
                        name: entry.name,
                        path: entryPath.join("/"),
                        isDirectory: entry.type === "directory",
                    };
                    if (entry.type === "file" && typeof entry.size === "number") {
                        file.size = entry.size;
                    }
                    if (entry.type === "directory" && Array.isArray(entry.children)) {
                        file.children = formatEntries(entry.children, entryPath);
                    }
                    return file;
                });
            };

            return formatEntries(data.entries, path);
        } else {
            console.error("Fetch error:", data.message);
            return [];
        }
    });

    return res;
};

export const fetchFiles = createAsyncThunk(
    "files/fetchFiles",
    async ({ uuid, path }: { uuid: string; path: string[] }) => {
      const files = await queryFiles(uuid, path);
      return files;
    }
  );

const initialState: File[] = [];

const fileSlice = createSlice({
    name: 'files',
    initialState,
    reducers: {},
    extraReducers: (builder) => {
        builder.addCase(fetchFiles.fulfilled, (state, action: PayloadAction<File[]>) => {
            return action.payload;
        });
    }
});

export const selectFiles = (state: { files: File[] }) => state.files;
export default fileSlice.reducer;
