import React from 'react';
import { ChordConfig } from '@/app/types/genomes';
import Box from '@mui/joy/Box';
import Slider from '@mui/joy/Slider';
import Checkbox from '@mui/joy/Checkbox';

interface ChordFormProps {
    onUpdate: (newConfig: ChordConfig) => void;
    defaultConfig: ChordConfig;
}

export const ChordForm: React.FC<ChordFormProps> = ({ onUpdate, defaultConfig }) => {
    const [config, setConfig] = React.useState(defaultConfig);

    const handleConfigChange = (key: keyof ChordConfig, value: any) => {
        const updatedConfig = { ...config, [key]: value };
        setConfig(updatedConfig);
        onUpdate(updatedConfig);
    };

    return (
        <div>
            <form>
            <Box 
                    sx={{
                        display: 'grid',
                        gridTemplateColumns: '20% 80%',
                        gap: 1,
                        alignItems: 'center',
                        width: '100%',
                        padding: 1,
                        fontSize: '0.65rem'
                    }}
                >
                    <label>Chord Padding</label>
                    <Slider 
                        valueLabelDisplay="auto" 
                        variant="solid" 
                        min={0} 
                        max={20} 
                        step={2} 
                        value={config.chordPadding} 
                        onChange={(e, value) => handleConfigChange('chordPadding', value)} 
                    />
                    <label>Opacity</label>
                    <Slider 
                        valueLabelDisplay="auto" 
                        variant="solid" 
                        min={0} 
                        max={1} 
                        step={0.1} 
                        value={config.opacity} 
                        onChange={(e, value) => handleConfigChange('opacity', value)} 
                    />
                    <label>Use Stroke</label>
                    <Checkbox 
                        checked={config.useStroke} 
                        onChange={(e) => handleConfigChange('useStroke', e.target.checked)}
                    />
                </Box>
            </form>
        </div>
    );
}