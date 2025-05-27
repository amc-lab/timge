"use client";
import { useState, useEffect } from "react";
import Segment from "./components/segment";
import Chords from "./components/chord";
import Bar from "./components/bar";
import { Track, TrackType } from "./config/track";
import Ring from "./components/ring";
import Line from "./components/line";
import Annotation from "./components/annotation";
import { AnnotationData as AnnotationType, GlobalConfig } from "../types/genomes";

interface TracksProps {
  tracks: Array<Track>;
  crossViewActionHandler?: any;
  id?: string;
  globalConfig?: GlobalConfig;
}

const Tracks = ({ tracks, crossViewActionHandler, id, globalConfig }: TracksProps) => {
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

  useEffect(() => {
    let minAvailableRadius = 160;

    const updatedTracks = [...tracks]
      .reverse()
      .map((track) => {
        if (track.trackType === TrackType.Karyotype) {
          const segmentInnerRadius = minAvailableRadius;
          const trackWidth = track.config.segmentTrackWidth;
          minAvailableRadius =
            segmentInnerRadius + trackWidth + track.config.segmentGridPadding;
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
    if (action === "generate_heatmap") {
      crossViewActionHandler(
        "propagate_dependencies",
        {
          viewId: id,
          dependencies: {
            segmentA: data.segmentA,
            segmentB: data.segmentB,
          }
        }
      );
    }
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
        // else if (track.trackType === TrackType.Highlight) {
        //   return (
        //     <Highlight
        //       key={index}
        //       divRef={track.data.divRef}
        //       segmentStartIdx={0}
        //       segmentEndIdx={0}
        //       segmentStartPos={2000}
        //       segmentEndPos={5000}
        //       segments={segmentData}
        //       globalConfig={track.data.globalConfig}
        //     />
        //   );
        // }
        return null;
      })}
    </>
  );
};

export default Tracks;
