import React from 'react';
import { ChordConfig } from '@/app/types/genomes';

interface ChordFormProps {
    onUpdate: (newConfig: ChordConfig) => void;
    defaultConfig: ChordConfig;
}

export const ChordForm: React.FC<ChordFormProps> = ({ onUpdate, defaultConfig }) => {
    return (
        <></>
    );
}