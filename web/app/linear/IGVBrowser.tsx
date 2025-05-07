'use client';
import igv from '../../node_modules/igv/dist/igv.esm.js';
import { useEffect, useRef } from 'react';
import { useAppSelector } from '@/store/hooks';

interface IGVBrowserProps {
    reference: string;
    trackFiles: string[];
}

const trackFormatMapping: { [key: string]: string } = {
    "bedgraph": "wig",
    "bigwig": "wig",
    "wig": "wig",
    "bp": "arc",
    "bed": "annotation",
};


const IGVBrowser = (props: IGVBrowserProps) => {
    const space = useAppSelector((state) => state.space);
    const HOST = process.env.NEXT_PUBLIC_DJANGO_HOST !== "http://127.0.0.1:8000" ? "https://timge.doc.ic.ac.uk/uploads/" + space.uuid + "/" : "http://localhost:3000/";
    const igvContainerRef = useRef<HTMLDivElement>(null);

    useEffect(() => {
        const igvDiv = igvContainerRef.current;

        if (!igvDiv) return;

        const options = {
            reference: {
                id: props.reference,
                name: props.reference,
                fastaURL: HOST + props.reference,
                indexURL: HOST + props.reference + ".fai",
                tracks: props.trackFiles.map((file) => {
                    const fileName = file.split("/").pop();
                    const fileType = fileName?.split(".").pop();
                    const trackType = trackFormatMapping[fileType || ""] || "annotation";
                    return {
                        name: fileName,
                        url: HOST + file,
                        type: trackType,
                        format: fileType,
                        displayMode: "EXPANDED",
                    };
                }),
            }
        };

        igv.createBrowser(igvDiv, options);

    }, [props.reference, props.trackFiles]);

    return (
        <div
            ref={igvContainerRef}
            style={{
                width: "100%",
            }}
        ></div>
    );
}

export default IGVBrowser;
