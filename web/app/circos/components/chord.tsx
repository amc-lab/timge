"use client";
import { useEffect } from "react";
import { Chord, ChordConfig, GlobalConfig } from "@/app/types/genomes";
import { useState } from "react";
import * as d3 from "d3";

interface ChordProps {
  data: {
    chords: Array<Chord>;
    globalConfig: GlobalConfig;
    divRef: any;
  };
  config: ChordConfig;
  segments: Array<any>;
  selectedSegments?: string[];
  idx: number;
  globalConfig?: GlobalConfig;
  dependencies?: any;
}

const Chords = ({ data, config, segments, selectedSegments, idx, globalConfig, dependencies }: ChordProps) => {
  const canvasRef = data.divRef;
  const chords = data.chords;

  useEffect(() => {
    console.log(globalConfig)
    if (!canvasRef.current) return;
    let svg = d3.select(canvasRef.current).select("svg");

    if (svg.empty()) {
      svg = d3
        .select(canvasRef.current)
        .append("svg")
        .attr("width", globalConfig.canvasWidth)
        .attr("height", globalConfig.canvasHeight);
    }

    const uniqueGroupClass = `group-${idx}`;

    svg.selectAll(`g.${uniqueGroupClass}`).remove();

    const group = svg
      .append("g")
      .attr("class", uniqueGroupClass)
      .attr(
        "transform",
        `translate(${globalConfig.canvasWidth / 2}, ${globalConfig.canvasHeight / 2})`,
      );

    const radius = config.outerRadius;
    const chord_padding = config.chordPadding;
    const chord_radius = radius - chord_padding;

    const minFilterScore = globalConfig.filterScore;
    // const maxFilterScore = config.maxFilterScore;
    let filteredChords = chords.filter((d) => d.score >= minFilterScore);
    filteredChords.sort((a, b) => b.score - a.score);

    console.log("Negative strand:", globalConfig.negativeStrand);
    let processedChords = filteredChords.map((chord) => {
      if (!globalConfig.negativeStrand) return chord;

      const sourceSegment = segments.find(s => s.chromosome === chord.source_chromosome);
      const targetSegment = segments.find(s => s.chromosome === chord.target_chromosome);
      if (!sourceSegment || !targetSegment) return chord;

      const len_x = sourceSegment.length;
      const len_y = targetSegment.length;

      return {
        ...chord,
        source_start: len_x - chord.source_end,
        source_end: len_x - chord.source_start,
        target_start: len_y - chord.target_end,
        target_end: len_y - chord.target_start,
      };
    });

    if (!dependencies) {
      dependencies = { chr: undefined, start: undefined, end: undefined };
    }
    const {chr, start, end} = dependencies;
    if (chr && chr !== "all" && start !== undefined && end !== undefined) {
      processedChords = processedChords.filter((d => {
        const sourceSegment = segments.find(segment => segment.chromosome === d.source_chromosome);
        const targetSegment = segments.find(segment => segment.chromosome === d.target_chromosome);
        if (!sourceSegment || !targetSegment) return false;
        const sourceStartAngle = sourceSegment.startAngle +
          (d.source_start / sourceSegment.length) *
          (sourceSegment.endAngle - sourceSegment.startAngle);
        const sourceEndAngle = sourceSegment.startAngle +
          (d.source_end / sourceSegment.length) *
          (sourceSegment.endAngle - sourceSegment.startAngle);
        const targetStartAngle = targetSegment.startAngle +
          (d.target_start / targetSegment.length) *
          (targetSegment.endAngle - targetSegment.startAngle);
        const targetEndAngle = targetSegment.startAngle +
          (d.target_end / targetSegment.length) *
          (targetSegment.endAngle - targetSegment.startAngle);
        return (
          (d.source_chromosome === chr && d.source_start >= start && d.source_end <= end) ||
          (d.target_chromosome === chr && d.target_start >= start && d.target_end <= end) ||
          (d.source_chromosome === chr && d.target_chromosome === chr &&
            ((d.source_start >= start && d.source_end <= end) ||
              (d.target_start >= start && d.target_end <= end)))
        );
      }));
    }

    const maxScore = d3.max(processedChords, (d) => d.score);
    const minScore = d3.min(processedChords, (d) => d.score);

    const colourScale = d3
      .scaleSequential()
      .domain([minScore, maxScore])
      .interpolator(d3.interpolateOrRd);

    const opacityScale = d3
      .scaleLinear()
      .domain([minScore, maxScore])
      .range([config.opacity, 1]);

    group
      .selectAll("path")
      .data(processedChords)
      .join("path")
      .attr(
        "d",
        d3
          .ribbon()
          .radius(chord_radius)
          .source((d) => ({
            startAngle:
              segments.find(
                (segment) => segment.chromosome === d.source_chromosome,
              )?.startAngle +
              (d.source_start /
                segments.find(
                  (segment) => segment.chromosome === d.source_chromosome,
                )?.length) *
                (segments.find(
                  (segment) => segment.chromosome === d.source_chromosome,
                )?.endAngle -
                  segments.find(
                    (segment) => segment.chromosome === d.source_chromosome,
                  )?.startAngle),
            endAngle:
              segments.find(
                (segment) => segment.chromosome === d.source_chromosome,
              )?.startAngle +
              (d.source_end /
                segments.find(
                  (segment) => segment.chromosome === d.source_chromosome,
                )?.length) *
                (segments.find(
                  (segment) => segment.chromosome === d.source_chromosome,
                )?.endAngle -
                  segments.find(
                    (segment) => segment.chromosome === d.source_chromosome,
                  )?.startAngle),
          }))
          .target((d) => ({
            startAngle:
              segments.find(
                (segment) => segment.chromosome === d.target_chromosome,
              )?.startAngle +
              (d.target_start /
                segments.find(
                  (segment) => segment.chromosome === d.target_chromosome,
                )?.length) *
                (segments.find(
                  (segment) => segment.chromosome === d.target_chromosome,
                )?.endAngle -
                  segments.find(
                    (segment) => segment.chromosome === d.target_chromosome,
                  )?.startAngle),
            endAngle:
              segments.find(
                (segment) => segment.chromosome === d.target_chromosome,
              )?.startAngle +
              (d.target_end /
                segments.find(
                  (segment) => segment.chromosome === d.target_chromosome,
                )?.length) *
                (segments.find(
                  (segment) => segment.chromosome === d.target_chromosome,
                )?.endAngle -
                  segments.find(
                    (segment) => segment.chromosome === d.target_chromosome,
                  )?.startAngle),
          })),
      )
      .attr("fill", (d) => {
        return colourScale(d.score);
      })
      .attr("opacity", (d) => {
        if (selectedSegments?.length > 0 && ! selectedSegments.includes(d.source_chromosome)) {
          return globalConfig.linkUnselectedOpacity || 0;
        }
        return globalConfig.linkSelectedOpacity || 0.6;
      })
      .attr("stroke", `${config.useStroke ? "black" : "none"}`)
      .on("mouseover", function () {
        d3.select(this).attr("opacity", 1);
      })
      .on("mouseout", function () {
        d3.select(this).attr("opacity", (d) => {
          if (selectedSegments?.length > 0 && ! selectedSegments.includes(d.source_chromosome)) {
            return globalConfig.linkUnselectedOpacity || 0;
          }
          return globalConfig.linkSelectedOpacity || 0.6;
        })
      });
  }, [canvasRef, chords, config, globalConfig, segments, selectedSegments, dependencies]);

  return null;
};

export default Chords;
