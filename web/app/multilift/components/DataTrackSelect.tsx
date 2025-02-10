import React from 'react';
import { Select, MenuItem } from '@mui/material';

interface DataTrackSelectProps {
    dataTracks: string[];
    selectedDataTrack: string;
    setSelectedDataTrack: (dataTrack: string) => void;
}

const DataTrackSelect = (props: DataTrackSelectProps) => {
    const { dataTracks, selectedDataTrack, setSelectedDataTrack } = props;
    return (
        <Select
            value={selectedDataTrack}
            onChange={(e) => setSelectedDataTrack(e.target.value)}
        >
            {dataTracks.map((dataTrack) => (
                <MenuItem key={dataTrack} value={dataTrack}>
                {dataTrack}
                </MenuItem>
            ))}
        </Select>
    );
}

export default DataTrackSelect;