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

export interface BarData {
  chromosome: string;
  start: number;
  end: number;
  value: number;
}

export interface RingData {
  chromosome: string;
  id: number;
  start: number;
  end: number;
}

export interface LineData {
  chromosome: string;
  chromStart: number;
  chromEnd: number;
  value: number;
}

interface AssemblyConfig {
  segmentPadding: number;
  axisLabelFontSize: number;
  showAxis: boolean;
  segmentInnerRadius: number;
  segmentTrackWidth: number;
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

interface BarConfig {
  innerRadius: number;
  trackWidth: number;
  trackPadding: number;
}

interface RingConfig {
  innerRadius: number;
  trackWidth: number;
  trackPadding: number;
  sequencePadding: number;
  tickLength: number;
  precision: number;
  metricPrefix: string;
  tickTextPadding: number;
  axisLabelFontSize: number;
  showAxis: boolean;
  gridPadding: number;
  hide: boolean;
}

interface LineConfig {
  innerRadius: number;
  trackWidth: number;
  trackPadding: number;
  colour?: string;
  hide: boolean;
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
