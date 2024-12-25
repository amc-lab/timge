export interface Assembly {
    chromosome: string;
    id: number;
    start: number;
    end: number;
    color?: string;
}

export interface Chord {
    source_chromosome: string;
    source_start: number;
    source_end: number;
    target_chromosome: string;
    target_start: number;
    target_end: number;
    color?: string;
}

interface AssemblyConfig {
    segmentPadding: number;
    axisLabelFontSize: number;
    showAxis: boolean;
    segmentInnerRadius: number;
    segmentOuterRadius: number;
    segmentGridPadding: number;
    canvasWidth: number;
    canvasHeight: number;
    tickLength: number;
    tickTextPadding: number;
    precision: number;
    useStroke: boolean;
}

interface ChordConfig {
    canvasWidth: number;
    canvasHeight: number;
}