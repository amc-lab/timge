import { useState } from 'react';
import Segment from './components/segment';
import Chords from './components/chord';
import { Track, TrackType } from './config/track';

interface TracksProps {
    tracks: Array<Track>;
}

const Tracks = (trackData: TracksProps) => {
    const [tracks, setTracks] = useState<Track[]>(trackData.tracks);
    const [segmentData, setSegmentData] = useState<any[]>([]);

    console.log(tracks);

    return (
        tracks.map((track, index) => {
            if (track.trackType === TrackType.karotype) {
                return <Segment key={index} data={track.data} config={track.config} onSegmentsCreated={(segData) => setSegmentData(segData)} />;
            }
            else if (track.trackType === TrackType.chord && segmentData.length > 0) {
                return <Chords key={index} data={track.data} config={track.config} segments={segmentData} />
            }
        })
    );
}

export default Tracks;