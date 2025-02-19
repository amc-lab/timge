import { useState, useEffect } from "react";
import Segment from "./components/segment";
import Chords from "./components/chord";
import Bar from "./components/bar";
import { Track, TrackType } from "./config/track";
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

  useEffect(() => {
    let minAvailableRadius = 160;

    const updatedTracks = [...tracks]
      .reverse()
      .map((track) => {
        if (track.trackType === TrackType.Karyotype) {
          const segmentInnerRadius = minAvailableRadius;
          const trackWidth = track.config.width;
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
        return null;
      })}
    </>
  );
};

export default Tracks;
