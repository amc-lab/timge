"use client";

import { useEffect, useRef, useState } from "react";
import * as d3 from "d3";
import { Assembly } from "../types/genomes";
import FileUpload from "../../components/FileUpload";
import Form from "../../components/Form";

const CanvasPage = () => {
    const canvasRef = useRef<HTMLDivElement>(null);
    const [segments, setSegments] = useState<Assembly[]>([]);
    const [axisLabelFontSize, setAxisLabelFontSize] = useState(10);
    const [showAxis, setShowAxis] = useState(true);
    const [segmentGridPadding, setSegmentGridPadding] = useState(3);
    const [tickLength, setTickLength] = useState(3);
    const [tickTextPadding, setTickTextPadding] = useState(3);
    const [segmentInnerRadius, setSegmentInnerRadius] = useState(200);
    const [segmentOuterRadius, setSegmentOuterRadius] = useState(250);
    const [canvasWidth, setCanvasWidth] = useState(600);
    const [canvasHeight, setCanvasHeight] = useState(600);
    const [segmentAnglePadding, setSegmentAnglePadding] = useState(0.02);
    const [precision, setPrecision] = useState(1);
    const [useStroke, setUseStroke] = useState(true);

    const colorPalette = d3.scaleSequential(d3.interpolateViridis);

    const handleFileUpload = (file: File) => {
        const reader = new FileReader();
        reader.onload = (event) => {
        try {
            const json = JSON.parse(event.target?.result as string);
            setSegments(json);
        } catch (error) {
            console.error("Invalid JSON file:", error);
            alert("Uploaded file is not valid JSON.");
        }
        };
        reader.readAsText(file);
    };

    useEffect(() => {
        if (canvasRef.current && segments.length > 0) {
            d3.select(canvasRef.current).select("svg").remove();

            const totalLength = segments.reduce(
                (acc, segment) => acc + (segment.end - segment.start),
                0
            );

            const svg = d3
                .select(canvasRef.current)
                .append("svg")
                .attr("viewBox", [-canvasWidth / 2, -canvasHeight / 2, canvasWidth, canvasHeight])
                .attr("width", canvasWidth)
                .attr("height", canvasHeight)
                .attr("style", "font: 10px sans-serif;");

            const chord = d3.chord()
                .padAngle(segmentAnglePadding)
                .sortSubgroups(d3.descending);

            const arc = d3.arc()
                .innerRadius(segmentInnerRadius)
                .outerRadius(segmentOuterRadius);

            const axisGridlineArc = d3.arc()
                .innerRadius(segmentOuterRadius + segmentGridPadding)
                .outerRadius(segmentOuterRadius + segmentGridPadding);

            const matrix = segments.map((segment, i) =>
                segments.map((_, j) => (i === j ? segment.end - segment.start : 0))
            );

            const chords = chord(matrix);

            const group = svg
                .append("g")
                .selectAll("g")
                .data(chords.groups)
                .join("g");

            group
                .append("path")
                .attr("fill", (d: { index: number; }) => colorPalette(0.3 + 0.6 * d.index / segments.length))
                .attr("d", arc)
                .attr("stroke", `${useStroke ? "black" : "none"}`)
                .on("mouseover", (event: any, d: { index: number; }) => {
                    d3.select(event.currentTarget).attr("filter", "brightness(0.95)");
                    }
                )
                .on("mouseout", (event: any, d: { index: number; }) => {
                    d3.select(event.currentTarget).attr("filter", "brightness(1)");
                    }
                )
                .append("title")

            group
                .append("text")
                .attr("x", 0)
                .attr("y", 0)
                .attr("transform", (d: { startAngle: number; endAngle: number; }) => {
                    const angle = (d.startAngle + d.endAngle) / 2;
                    const x = Math.sin(angle) * (segmentInnerRadius + (segmentOuterRadius - segmentInnerRadius) / 2);
                    const y = -Math.cos(angle) * (segmentInnerRadius + (segmentOuterRadius - segmentInnerRadius) / 2);
                    return `translate(${x}, ${y})`;
                })
                .attr("text-anchor", "middle")
                .attr("alignment-baseline", "middle")
                .attr("font-size", `10`)
                .attr("fill", "white")
                .text((d: { index: number; }) => segments[d.index].chromosome);

            const precision_vals = [100, 250, 500]
            const standardIntervalSizes = [500, 1000, 2500, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000]
            const expectedTickPrecision = totalLength / precision_vals[precision];
            const minorTickInterval = d3.bisect(standardIntervalSizes, expectedTickPrecision) === 0 ? standardIntervalSizes[0] : standardIntervalSizes[d3.bisect(standardIntervalSizes, expectedTickPrecision) - 1];
            const tickInterval = minorTickInterval * 5;

            if (showAxis) {
                group
                    .append("path")
                    .attr("d", axisGridlineArc)
                    .attr("stroke", "black")

                group
                    .append("g")
                    .selectAll("g")
                    .data((d: { index: string | number; startAngle: any; endAngle: any; }) => {
                    const segment = segments[d.index];
                    const scale = d3.scaleLinear()
                        .domain([segment.start, segment.end])
                        .range([d.startAngle, d.endAngle]);

                    const majorTicks = d3.range(segment.start, segment.end, tickInterval).map((value: number) => ({
                        value: value - segment.start == 0 ? null : d3.formatPrefix(",.0", 1e3)(value),
                        angle: scale(value),
                        isMajor: true,
                    }));

                    const minorTicks = d3.range(segment.start, segment.end, minorTickInterval).filter((value: number) => value % tickInterval !== 0).map((value: any) => ({
                        value: null,
                        angle: scale(value),
                        isMajor: false,
                    }));

                    return [...majorTicks, ...minorTicks];
                    })
                    .join("g")
                    .attr("transform", (d: { angle: number; }) => `rotate(${(d.angle * 180) / Math.PI - 90}) translate(${segmentOuterRadius + segmentGridPadding}, 0)`)
                    .call((g) =>
                    g
                        .append("line")
                        .attr("stroke", "currentColor")
                        .attr("x2", (d: { isMajor: any; }) => (d.isMajor ? tickLength * 2 : tickLength))
                    )
                    .call((g) =>
                    g
                        .filter((d: { isMajor: any; }) => d.isMajor)
                        .append("text")
                        .attr("x", 2 * tickLength + tickTextPadding)
                        .attr("dy", "0.35em")
                        .attr("font-size", `${axisLabelFontSize}`)
                        .attr("transform", (d: { angle: number; }) =>
                        d.angle > Math.PI ? `rotate(180) translate(-${4 * tickLength + 2 * tickTextPadding})` : null
                        )
                        .attr("text-anchor", (d: { angle: number; }) => (d.angle > Math.PI ? "end" : null))
                        .text((d: { value: { toLocaleString: () => any; }; }) => d.value?.toLocaleString())
                    );
            }
        }
    }, [segments, segmentInnerRadius, segmentOuterRadius, segmentGridPadding, showAxis, axisLabelFontSize, tickLength, tickTextPadding, canvasWidth, canvasHeight, segmentAnglePadding, precision, useStroke]);

    return (
        <div>
        <h1>D3 Chord Diagram</h1>
        <FileUpload onFileUpload={handleFileUpload} />
        <div className="svg-container" style={{display: "flex", justifyContent: "center", alignItems: "center"}}>
        <div ref={canvasRef} style={{ border: "1px solid black" }}></div>
        </div>
        {segments.length > 0 && 
            <Form onUpdate={
                (
                    segmentPadding, 
                    axisLabelFontSize, 
                    showAxis, 
                    segmentInnerRadius, 
                    segmentOuterRadius, 
                    segmentGridPadding, 
                    canvasWidth,
                    canvasHeight,
                    tickLength,
                    tickTextPadding,
                    precision,
                    useStroke
                ) => {
                setSegmentAnglePadding(segmentPadding);
                setAxisLabelFontSize(axisLabelFontSize);
                setShowAxis(showAxis);
                setSegmentInnerRadius(segmentInnerRadius);
                setSegmentOuterRadius(segmentOuterRadius);
                setSegmentGridPadding(segmentGridPadding);
                setCanvasWidth(canvasWidth);
                setCanvasHeight(canvasHeight);
                setTickLength(tickLength);
                setTickTextPadding(tickTextPadding);
                setPrecision(precision);
                setUseStroke(useStroke);
            }} />
        }
        </div>
    );
};

export default CanvasPage;
