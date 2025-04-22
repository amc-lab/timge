import React from 'react';
import { Select, ListItem, Option } from '@mui/joy';
import { Card } from '@mui/joy';
import { Typography } from '@mui/joy';

interface DataTrackSelectProps {
    dataTracks: File[];
    genome: string;
    setSelectedDataTrack: (dataTrack: File, genome: string) => void;
}

const DataTrackSelect = (props: DataTrackSelectProps) => {
    const { dataTracks, genome, setSelectedDataTrack } = props;
    const [selectedValues, setSelectedValues] = React.useState<File[]>([]);

    return (
        <Card
            sx={{
                width: "calc(50% - 0.5em)",
                display: "flex",
                marginTop: "0.5em",
            }}
        >
        <Select
            defaultValue={[]}
            placeholder={`Select track files for ${genome}`} 
            sx={{
                width: "100%",
                display: "flex",
                marginTop: "0.5em",
                marginBottom: "0.5em",
            }}
            multiple
            onChange={(event, value) => {
                for (const dataTrack of selectedValues) {
                    if (!value.includes(dataTrack)) {
                        setSelectedDataTrack(dataTrack, null);
                    }
                }
                for (const dataTrack of value) {
                    setSelectedDataTrack(dataTrack, genome);
                }
                setSelectedValues(value);
            }}
            slotProps={{
                listbox: {
                  sx: {
                    zIndex: 1400,
                  },
                },
              }}
        >
            {dataTracks.map((dataTrack) => (
                <Option key={dataTrack.name} value={dataTrack}>
                    <Typography>{dataTrack.name}</Typography>
                </Option>
            ))}
        </Select>
        </Card>
    );
}

export default DataTrackSelect;