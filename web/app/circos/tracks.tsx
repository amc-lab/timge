import { useState } from 'react';
import Segment from './components/segment';
import Chords from './components/chord';
import { Track, TrackType } from './config/track';

interface TracksProps {
    tracks: Array<Track>;
}

const Tracks = (trackData: TracksProps) => {
    const [segmentData, setSegmentData] = useState<any[]>([]);

    const onSegmentsCreated = (segData: any[]) => {
        setSegmentData(segData);
    }

    return (
        trackData.tracks.map((track, index) => {
            if (track.trackType === TrackType.Karyotype) {
                return <Segment key={index} data={track.data} config={track.config} onSegmentsCreated={onSegmentsCreated} />;
            }
            else if (track.trackType === TrackType.Chord && segmentData.length > 0) {
                return <Chords key={index} data={track.data} config={track.config} segments={segmentData} />
            }
        })
    );
}

export default Tracks;