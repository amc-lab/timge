  export const uploadTrackFiles = async (space, files: File[]) => {
    const host = process.env.NEXT_PUBLIC_DJANGO_HOST;
    const formData = new FormData();

    files.forEach((track) => {
      formData.append("track_files", track);
    });
    formData.append("uuid", space.uuid);
    formData.append("path", space.config.working_directory);


    const response = await fetch(`${host}/api/timge/upload_tracks/`, {
      method: "POST",
      body: formData,
    });

    const data = await response.json();

    if (response.ok) {
      console.log("Files uploaded successfully", data);
      return data;
    } else {
      throw new Error(data.message || "Upload failed");
    }
  };

  const fileFormatMapping = {
    "fasta": "karyotype",
    "bed": "line",
    "bedgraph": "line",
    "bedpe": "link",
    "fa": "karyotype",
  }

  export const getTrackFiles = async (space: any, useWorkingDirectory: boolean): Promise<any[]> => {
    try {
      const host = process.env.NEXT_PUBLIC_DJANGO_HOST;
      const response = await fetch(`${host}/api/timge/get_files_in_path/`, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({
          uuid: space.uuid,
          path: useWorkingDirectory ? space.config.working_directory : [],
        }),
      });
  
      const data = await response.json();
  
      if (data.status === "success") {
        return data.entries.map((file) => ({
          name: file.name,
          type: file.type,
          size: file.size,
        }));
      } else {
        console.error("Error fetching tracks:", data.message);
        return [];
      }
    } catch (err) {
      console.error("Request failed:", err);
      return [];
    }
  };
  