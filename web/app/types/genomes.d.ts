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
    tickLength: number;
    tickTextPadding: number;
    precision: number;
    useStroke: boolean;
    metricPrefix: string;
}

interface ChordConfig {
    chordPadding: number;
    opacity: number;
    colour?: string;
    useStroke: boolean;
    outerRadius: number;
}

interface GlobalConfig {
    canvasWidth: number;
    canvasHeight: number;
}

interface SegmentDetails {
    chromosome: string;
    start: number;
    end: number;
    length: number;
    innerRadius: number;
}
