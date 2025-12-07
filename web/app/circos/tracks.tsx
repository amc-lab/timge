"use client";
import { useState, useEffect } from "react";
import Segment from "./components/segment";
import Chords from "./components/chord";
import Bar from "./components/bar";
import { Track, TrackType } from "./config/track";
import Ring from "./components/ring";
import Line from "./components/line";
import Annotation from "./components/annotation";
import Highlight from "./components/highlight";
import { AnnotationData as AnnotationType, Chord, GlobalConfig } from "../types/genomes";
import { publishCrossViewEvent } from "@/app/utils/crossViewEvents";

interface TracksProps {
  tracks: Array<Track>;
  id?: string;
  globalConfig?: GlobalConfig;
  dependencies?: any;
}

const DEFAULT_HEATMAP_RESOLUTION = 100;

const Tracks = ({ tracks, id, globalConfig, dependencies }: TracksProps) => {
  const [segmentData, setSegmentData] = useState<any[]>([]);
  const [selectedSegments, setSelectedSegments] = useState<string[]>([]);
  const [trackData, setTrackData] = useState<Array<Track>>([]);

  const onSegmentsCreated = (segData: any[]) => {
    setSegmentData(segData);
  };

  useEffect(() => {
    setTrackData(tracks);
  }, [tracks]);

  const [totalRadius, setTotalRadius] = useState(160);

  useEffect(() => {
    console.log(totalRadius);
  }
  , [totalRadius]);

  const setHighlightConfig = (tracks: Array<Track>, newHighlightConfig: any) => {
    const highlightTrack = tracks.find(track => track.trackType === TrackType.Highlight);
    if (highlightTrack) {
      highlightTrack.config = {
        ...highlightTrack.config,
        ...newHighlightConfig
      };
    }
  };

  const generateRNAfold = (d: Chord) => {
    const reference =  tracks.find(track => track.trackType === TrackType.Karyotype)?.name || "";
    if (!reference) {
      console.error("No karyotype track found to generate RNAfold.");
      return;
    }
    const segmentA = d.source_chromosome;
    const segmentB = d.target_chromosome;
    const segmentAStart = d.source_start;
    const segmentAEnd = d.source_end;
    const segmentBStart = d.target_start;
    const segmentBEnd = d.target_end;
    
    publishCrossViewEvent("GENERATE_RNAFOLD", {
      reference,
      segmentA: d.source_chromosome,
      segmentB: d.target_chromosome,
      segmentAStart: d.source_start,
      segmentAEnd: d.source_end,
      segmentBStart: d.target_start,
      segmentBEnd: d.target_end,
    });
  }

  useEffect(() => {
    let minAvailableRadius = 160;
    let karyotypeWidth = 0;

    const updatedTracks = [...tracks]
      .reverse()
      .map((track) => {
        if (track.trackType === TrackType.Karyotype) {
          const segmentInnerRadius = minAvailableRadius;
          const trackWidth = track.config.segmentTrackWidth;
          minAvailableRadius =
            segmentInnerRadius + trackWidth + track.config.segmentGridPadding;
          setHighlightConfig(tracks, {
            innerRadius: segmentInnerRadius - 5,
            width: trackWidth + 10,
          });
          return {
            ...track,
            config: {
              ...track.config,
              segmentInnerRadius,
            },
          };
        } else if (track.trackType === TrackType.Bar) {
          const barInnerRadius = minAvailableRadius;
          minAvailableRadius =
            barInnerRadius +
            track.config.trackWidth +
            track.config.trackPadding;
          return {
            ...track,
            config: {
              ...track.config,
              innerRadius: barInnerRadius,
            },
          };
        } else if (track.trackType === TrackType.Chord) {
          minAvailableRadius = track.config.outerRadius;
          return track;
        } else if (track.trackType === TrackType.Ring) {
          if (track.config.hide) {
            return track;
          }
          const ringInnerRadius = minAvailableRadius;
          minAvailableRadius =
            ringInnerRadius +
            track.config.trackWidth +
            track.config.trackPadding;
          if (track.config.showAxis) {
            minAvailableRadius += track.config.tickLength +
             track.config.tickTextPadding +
              track.config.axisLabelFontSize * 3 +
              track.config.gridPadding;
          }
          return {
            ...track,
            config: {
              ...track.config,
              innerRadius: ringInnerRadius,
            },
          };
        } else if (track.trackType === TrackType.Line) {
          const lineInnerRadius = minAvailableRadius;
          minAvailableRadius +=
            track.config.trackWidth +
            track.config.trackPadding;
          return {
            ...track,
            config: {
              ...track.config,
              innerRadius: lineInnerRadius,
            },
          };
        }
        else if (track.trackType === TrackType.Annotation) {
          const annotationInnerRadius = minAvailableRadius;
          console.log(findMaxOverlaps(track.data.annotations));
          minAvailableRadius =
            annotationInnerRadius +
            findMaxOverlaps(track.data.annotations) * (
            track.config.trackWidth +
            track.config.trackPadding +
            track.config.textFontSize + 
            track.config.textPadding);
          return {
            ...track,
            config: {
              ...track.config,
              innerRadius: annotationInnerRadius,
            },
          };
        }
        setTotalRadius(minAvailableRadius);
        return track;
      })
      .reverse();

    setTrackData(updatedTracks);
  }, [tracks]);

    const findMaxOverlaps = (annotations: AnnotationType[]): number => {
      const intervals = [];
      annotations.forEach((annotation) => {
        intervals.push({ start: annotation.chromStart, end: annotation.chromEnd });
      });
      intervals.sort((a, b) => a.start);
      let maxOverlaps = 1;
      let currentOverlaps = 1;
      let prevEnd = intervals[0].end;

      for (let i = 1; i < intervals.length; i++) {
        if (intervals[i].start < prevEnd) {
          currentOverlaps++;
          maxOverlaps = Math.max(maxOverlaps, currentOverlaps);
        } else {
          currentOverlaps = 1;
        }
        prevEnd = Math.max(prevEnd, intervals[i].end);
      }
      return maxOverlaps;
    };

  const onSelectSegments = (selectedSegments: string[]) => {
    console.log("Selected segments:", selectedSegments);
    setSelectedSegments(selectedSegments);
  }

  const onCustomAction = (action: string, data: any) => {
    console.log("Custom action triggered:", action, data);
    if (action !== "generate_heatmap") return;

    const segments = (data?.segments ?? []) as Array<{ id: string; length: number }>;
    if (segments.length < 2) {
      console.warn("Generate heatmap requires exactly two segments.");
      return;
    }

    const [segmentA, segmentB] = segments;

    if (!segmentA || !segmentB) {
      console.warn("Generate heatmap invoked without both segments.");
      return;
    }

    if (id) {
      publishCrossViewEvent("PROPAGATE_DEPENDENCIES", {
        viewId: id,
        dependencies: { segmentA, segmentB },
      });
    }

    const referenceTrack = trackData.find((track) => track.trackType === TrackType.Karyotype)
      ?.name;
    const contactTrack = trackData.find((track) => track.trackType === TrackType.Chord)?.name;

    if (!referenceTrack || !contactTrack) {
      console.warn(
        "Unable to generate heatmap from circos view without both reference and interaction tracks.",
      );
      return;
    }

    publishCrossViewEvent("GENERATE_HEATMAP", {
      reference: referenceTrack,
      track: contactTrack,
      segmentA,
      segmentB,
      resolution: DEFAULT_HEATMAP_RESOLUTION,
    });
  };

  return (
    <>
      {trackData.map((track, index) => {
        if (track.trackType === TrackType.Karyotype) {
          return (
            <>
            <Segment
              key={index}
              data={track.data}
              config={track.config}
              onSegmentsCreated={onSegmentsCreated}
              onSelectSegments={onSelectSegments}
              onCustomAction={onCustomAction}
              idx={id + index}
            />
            {/* {
              segmentData.length > 0 && (
                <Highlight
                  divRef={track.data.divRef}
                  key={index + "highlight"}
                  segmentStartIdx={0}
                  segmentEndIdx={0}
                  segmentStartPos={0}
                  segmentEndPos={500}
                  segments={segmentData}
                  globalConfig={globalConfig}
                  innerRadius={track.config.segmentInnerRadius}
                  width={track.config.segmentTrackWidth}
            />
              )
            } */}
            </>
          );
        } else if (track.trackType === TrackType.Bar) {
          return (
            <Bar
              key={index}
              data={track.data}
              config={track.config}
              segments={segmentData}
              idx={index}
            />
          );
        } else if (
          track.trackType === TrackType.Chord &&
          segmentData.length > 0
        ) {
          return (
            <Chords
              key={index}
              data={track.data}
              config={track.config}
              segments={segmentData}
              selectedSegments={selectedSegments}
              idx={index}
              globalConfig={globalConfig}
              dependencies={dependencies}
              onChordClick={(d) => {generateRNAfold(d)}}
            />
          );
        }
        else if (track.trackType === TrackType.Ring) {
          return (
            <Ring
              key={index}
              data={track.data}
              config={track.config}
              segments={segmentData}
              idx={index}
            />
          )
        }
        else if (track.trackType === TrackType.Line) {
          return (
            <Line
              key={index}
              data={track.data}
              config={track.config}
              segments={segmentData}
              idx={index}
              trackName={track.name}
            />
          )
        }
        else if (track.trackType === TrackType.Annotation) {
          return (
            <Annotation
              key={index}
              data={track.data}
              config={track.config}
              segments={segmentData}
              idx={index}
            />
          );
        }
        else if (track.trackType === TrackType.Highlight) {
          return (
            <Highlight
              key={index}
              segments={segmentData}
              config={track.config}
              data={track.data}
              dependencies={dependencies}
              globalConfig={globalConfig}
            />
          );
        }
        return null;
      })}
    </>
  );
};

export default Tracks;
