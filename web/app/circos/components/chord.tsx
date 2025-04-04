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
  selectedSegments?: Set<string>;
  idx: number;
}

const Chords = ({ data, config, segments, selectedSegments, idx }: ChordProps) => {
  const canvasRef = data.divRef;
  const chords = data.chords;
  const globalConfig = data.globalConfig;

  useEffect(() => {
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

    const minAcceptedScore = 100;
    const filteredChords = chords.filter((d) => d.score > minAcceptedScore);

    // get maximum score from all the segments
    const maxScore = d3.max(filteredChords, (d) => d.score);
    const minScore = d3.min(filteredChords, (d) => d.score);

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
      .data(filteredChords)
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
        if (selectedSegments?.size > 0 && ! selectedSegments.has(d.source_chromosome)) {
          return 0.2;
        }
        return config.opacity;
      })
      .attr("stroke", `${config.useStroke ? "black" : "none"}`)
      .on("mouseover", function () {
        d3.select(this).attr("opacity", 1);
      })
      .on("mouseout", function () {
        d3.select(this).attr("opacity", (d) => {
          if (selectedSegments?.size > 0 && ! selectedSegments.has(d.source_chromosome)) {
            return 0.2;
          }
          return config.opacity;
        })
      });
  }, [canvasRef, chords, config, globalConfig, segments, selectedSegments]);

  return null;
};

export default Chords;
