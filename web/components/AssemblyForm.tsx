import React, { useState } from 'react';
import Checkbox from '@mui/joy/Checkbox';
import Box from '@mui/joy/Box';
import Slider from '@mui/joy/Slider';
import Card from '@mui/joy/Card';
import { AssemblyConfig } from '@/app/types/genomes';

interface AssemblyFormProps {
    onUpdate: (newConfig: AssemblyConfig) => void;
    defaultConfig: AssemblyConfig;
}

export const AssemblyForm: React.FC<AssemblyFormProps> = ({ onUpdate, defaultConfig }) => {
    const [config, setConfig] = useState(defaultConfig);

    const handleConfigChange = (key: keyof AssemblyConfig, value: any) => {
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
                    <label>Segment Padding</label>
                    <Slider 
                        valueLabelDisplay="auto" 
                        variant="solid" 
                        min={0} 
                        max={0.1} 
                        step={0.01} 
                        value={config.segmentPadding} 
                        onChange={(e, value) => handleConfigChange('segmentPadding', value)} 
                    />
                    <label>Axis Font Size</label>
                    <Slider 
                        valueLabelDisplay="auto" 
                        variant="solid" 
                        min={6} 
                        max={16} 
                        step={1} 
                        value={config.axisLabelFontSize} 
                        onChange={(e, value) => handleConfigChange('axisLabelFontSize', value)} 
                    />
                    <label>Show Axis</label>
                    <Checkbox 
                        label="" 
                        checked={config.showAxis} 
                        onChange={(e) => handleConfigChange('showAxis', e.target.checked)} 
                    />
                    <label>Inner Radius</label>
                    <Slider 
                        valueLabelDisplay="auto" 
                        variant="solid" 
                        min={0} 
                        max={500} 
                        step={10} 
                        value={config.segmentInnerRadius} 
                        onChange={(e, value) => handleConfigChange('segmentInnerRadius', value)} 
                    />
                    <label>Outer Radius</label>
                    <Slider 
                        valueLabelDisplay="auto" 
                        variant="solid" 
                        min={0} 
                        max={500} 
                        step={10} 
                        value={config.segmentOuterRadius} 
                        onChange={(e, value) => handleConfigChange('segmentOuterRadius', value)} 
                    />
                    <label>Grid Padding</label>
                    <Slider 
                        valueLabelDisplay="auto" 
                        variant="solid" 
                        min={1} 
                        max={5} 
                        step={1} 
                        value={config.segmentGridPadding} 
                        onChange={(e, value) => handleConfigChange('segmentGridPadding', value)} 
                    />
                    <label>Canvas Width</label>
                    <Slider 
                        valueLabelDisplay="auto" 
                        variant="solid" 
                        min={0} 
                        max={1000} 
                        step={10} 
                        value={config.canvasWidth} 
                        onChange={(e, value) => handleConfigChange('canvasWidth', value)} 
                    />
                    <label>Canvas Height</label>
                    <Slider 
                        valueLabelDisplay="auto" 
                        variant="solid" 
                        min={0} 
                        max={1000} 
                        step={10} 
                        value={config.canvasHeight} 
                        onChange={(e, value) => handleConfigChange('canvasHeight', value)} 
                    />
                    <label>Tick Length</label>
                    <Slider 
                        valueLabelDisplay="auto" 
                        variant="solid" 
                        min={0} 
                        max={10} 
                        step={1} 
                        value={config.tickLength} 
                        onChange={(e, value) => handleConfigChange('tickLength', value)} 
                    />
                    <label>Text Padding</label>
                    <Slider 
                        valueLabelDisplay="auto" 
                        variant="solid" 
                        min={0} 
                        max={10} 
                        step={1} 
                        value={config.tickTextPadding} 
                        onChange={(e, value) => handleConfigChange('tickTextPadding', value)} 
                    />
                    <label>Precision</label>
                    <Slider 
                        valueLabelDisplay="auto" 
                        variant="solid" 
                        min={0} 
                        max={2} 
                        step={1} 
                        value={config.precision} 
                        onChange={(e, value) => handleConfigChange('precision', value)} 
                    />
                    <label>Use Stroke</label>
                    <Checkbox 
                        label="" 
                        checked={config.useStroke} 
                        onChange={(e) => handleConfigChange('useStroke', e.target.checked)} 
                    />
                </Box>
            </form>
        </div>
    );
}