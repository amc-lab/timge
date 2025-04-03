"use client";
import { useState, useEffect } from "react";
import Segment from "./components/segment";
import Chords from "./components/chord";
import Bar from "./components/bar";
import { Track, TrackType } from "./config/track";
import Highlight from "./components/highlight"; // Ensure this is the correct path to your Highlight component
import Ring from "./components/ring";
import Line from "./components/line";

interface TracksProps {
  tracks: Array<Track>;
}

const Tracks = ({ tracks }: TracksProps) => {
  const [segmentData, setSegmentData] = useState<any[]>([]);
  const [trackData, setTrackData] = useState<Array<Track>>([]);

  const onSegmentsCreated = (segData: any[]) => {
    setSegmentData(segData);
  };

  useEffect(() => {
    console.log("Setting track data", tracks);
    setTrackData(tracks);
  }, [tracks]);

  const [totalRadius, setTotalRadius] = useState(160);

  useEffect(() => {
    console.log(totalRadius);
  }
  , [totalRadius]);

  useEffect(() => {
    console.log("DEBUG: Setting track data", tracks);
    let minAvailableRadius = 160;

    const updatedTracks = [...tracks]
      .reverse()
      .map((track) => {
        if (track.trackType === TrackType.Karyotype) {
          const segmentInnerRadius = minAvailableRadius;
          const trackWidth = track.config.segmentTrackWidth;
          minAvailableRadius =
            segmentInnerRadius + trackWidth + track.config.segmentGridPadding;
          console.log("DEBUG: Karyotype track width", trackWidth);
          console.log("DEBUG: Karyotype segment inner radius", segmentInnerRadius);
          console.log("DEBUG: Karyotype segment grid padding", track.config.segmentGridPadding);
          console.log("DEBUG: Karyotype minimum available radius", minAvailableRadius);
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
          console.log("DEBUG: Bar minimum available radius", minAvailableRadius);
          return {
            ...track,
            config: {
              ...track.config,
              innerRadius: barInnerRadius,
            },
          };
        } else if (track.trackType === TrackType.Chord) {
          minAvailableRadius = track.config.outerRadius;
          console.log("DEBUG: Chord minimum available radius", minAvailableRadius);
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
          console.log("DEBUG: Line inner radius", lineInnerRadius);
          console.log("DEBUG: Line track width", track.config.trackWidth);
          console.log("DEBUG: Line track padding", track.config.trackPadding);
          minAvailableRadius +=
            track.config.trackWidth +
            track.config.trackPadding;
          console.log("DEBUG: Line minimum available radius", minAvailableRadius);
          return {
            ...track,
            config: {
              ...track.config,
              innerRadius: lineInnerRadius,
            },
          };
        }
        setTotalRadius(minAvailableRadius);
        return track;
      })
      .reverse();

    setTrackData(updatedTracks);
  }, [tracks]);

  return (
    <>
      {trackData.map((track, index) => {
        if (track.trackType === TrackType.Karyotype) {
          return (
            <Segment
              key={index}
              data={track.data}
              config={track.config}
              onSegmentsCreated={onSegmentsCreated}
              idx={index}
            />
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
              idx={index}
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
            />
          )
        }
        else if (track.trackType === TrackType.Highlight) {
          return (
            <Highlight
              key={index}
              divRef={track.data.divRef}
              segmentStartIdx={0}
              segmentEndIdx={0}
              segmentStartPos={2000}
              segmentEndPos={5000}
              totalRadius={totalRadius}
              segments={segmentData}
              globalConfig={track.data.globalConfig}
            />
          );
        }
        return null;
        return null;
      })}
    </>
  );
};

export default Tracks;
