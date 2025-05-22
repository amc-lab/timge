"use client";

import { useEffect, useRef } from "react";
import { useAppSelector } from "@/store/hooks";

declare global {
    interface Window {
      igv: any;
    }
  }

const IGVBrowser = ({ reference, trackFiles, viewConfig }) => {
  const space = useAppSelector((state) => state.space);
  const igvContainerRef = useRef(null);

  const HOST =
    process.env.NEXT_PUBLIC_DJANGO_HOST !== "http://127.0.0.1:8000"
      ? "https://timge.doc.ic.ac.uk/uploads/" + space.uuid + "/"
      : "http://localhost:3000/";

  useEffect(() => {
    const loadScript = () => {
      return new Promise<void>((resolve, reject) => {
        if (window.igv) return resolve();

        const script = document.createElement("script");
        script.src = "https://cdn.jsdelivr.net/npm/igv@3.0.2/dist/igv.min.js";
        script.async = true;
        script.onload = () => resolve();
        script.onerror = () => reject("IGV.js failed to load");
        document.body.appendChild(script);
      });
    };

    const initializeBrowser = async () => {
      if (!igvContainerRef.current) return;

      await loadScript();

      const options = {
        reference: {
          id: reference,
          name: reference,
          fastaURL: HOST + reference,
          indexURL: HOST + reference + ".fai",
          tracks: trackFiles.map((file) => {
            const fileName = file.split("/").pop();
            const fileType = fileName?.split(".").pop();
            const trackType =
              { bedgraph: "wig", bigwig: "wig", wig: "wig", bp: "arc", bed: "annotation" }[
                fileType || ""
              ] || "annotation";
            return {
              name: fileName,
              url: HOST + file,
              type: trackType,
              format: fileType,
            };
          }),
        },
      };

      // @ts-ignore
      const browser = await window.igv.createBrowser(igvContainerRef.current, options);
    };

    initializeBrowser();
  }, []);

  return <div ref={igvContainerRef} style={{ width: "100%" }} />;
};

export default IGVBrowser;