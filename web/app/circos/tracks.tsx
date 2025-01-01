import { useState, useEffect } from 'react';
import Segment from './components/segment';
import Chords from './components/chord';
import { Track, TrackType } from './config/track';

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
        setTrackData(tracks);
    }, [tracks]);

    useEffect(() => {
        let minAvailableRadius = 160;

        const updatedTracks = [...tracks]
            .reverse()
            .map((track) => {
                if (track.trackType === TrackType.Karyotype) {
                    const segmentInnerRadius = minAvailableRadius;
                    const trackWidth = track.config.width || 50;
                    minAvailableRadius = segmentInnerRadius + trackWidth + (track.config.segmentGridPadding || 0);
                    return {
                        ...track,
                        config: {
                            ...track.config,
                            segmentInnerRadius,
                        },
                    };
                } else {
                    minAvailableRadius = track.config.outerRadius || minAvailableRadius;
                    return track;
                }
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
                        />
                    );
                } else if (track.trackType === TrackType.Chord && segmentData.length > 0) {
                    return (
                        <Chords
                            key={index}
                            data={track.data}
                            config={track.config}
                            segments={segmentData}
                        />
                    );
                }
                return null;
            })}
        </>
    );
};

export default Tracks;
