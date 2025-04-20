  export const uploadTrackFiles = async(space, files: File[]) => {
    const host = process.env.NEXT_PUBLIC_DJANGO_HOST;
    let formData = new FormData();
    
    files.forEach((track) => {
      formData.append("track_files", track);
    });
    formData.append("uuid", space.uuid);
    fetch(`${host}/api/timge/upload_tracks/`, {
      method: "POST",
      body: formData,
    })
    .then((response) => response.json())
    .then((data) => {
      console.log("Files uploaded successfully", data);
      const trackFiles = files.map((file) => ({
        name: file.name,
        data: file,
        type: file.type,
        size: file.size,
        trackType: fileFormatMapping[file.name.split('.').pop()],
      }));
    })
    .catch((error) => {
      console.error("Error uploading files", error);
    });
  }

  const fileFormatMapping = {
    "fasta": "karyotype",
    "bed": "line",
    "bedgraph": "line",
    "bedpe": "link",
    "fa": "karyotype",
  }

  export const getTrackFiles = (space) => {
    const host = process.env.NEXT_PUBLIC_DJANGO_HOST;
    fetch(`${host}/api/timge/get_files_in_path/?uuid=${space.uuid}`, {
      method: "GET",
    })
      .then((response) => response.json())
      .then((data) => {
        if (data.status === "success") {
          const trackFiles = data.files.map((file) => ({
            name: file.name,
            size: file.size,
          }));
          return trackFiles;
        } else {
          console.error("Error fetching tracks:", data.message);
        }
      })
      .catch((err) => {
        console.error("Request failed:", err);
      });
  };