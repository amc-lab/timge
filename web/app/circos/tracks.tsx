import { useState } from 'react';
import { Track, TrackType } from '../types/genomes';
import Segment from './components/segment';
import Chords from './components/chord';

const Tracks = () => {
    const [tracks, setTracks] = useState(Array<Track>());

    return (
        tracks.map((track, index) => {
            if (track.trackType === TrackType.karotype) {
                return <Segment key={index} data={track.data} config={track.config} />;
            }
            else if (track.trackType === TrackType.chord) {
                return <Chords key={index} data={track.data} config={track.config} />;
            }
        })
    );
}