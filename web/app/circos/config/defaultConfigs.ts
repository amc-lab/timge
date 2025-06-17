import { text } from "stream/consumers";

export const defaultAssemblyConfig = {
  segmentPadding: 0.02,
  axisLabelFontSize: 10,
  showAxis: true,
  segmentInnerRadius: 200,
  segmentTrackWidth: 40,
  segmentGridPadding: 0,
  tickLength: 5,
  tickTextPadding: 2,
  precision: 1,
  useStroke: true,
  metricPrefix: "k",
};

export const defaultChordConfig = {
  chordPadding: 10,
  opacity: 0.6,
  color: "blue",
  useStroke: true,
  outerRadius: defaultAssemblyConfig.segmentInnerRadius,
  minFilterScore: 0,
  maxFilterScore: 1000,
};

export const defaultHighlightConfig = {
  innerRadius: 195,
  width: 50,
  showHighlight: false,
};

export const defaultBarConfig = {
  innerRadius: 200,
  trackPadding: 10,
  trackWidth: 24,
};

export const defaultRingConfig = {
  innerRadius: 200,
  trackPadding: 10,
  trackWidth: 24,
  axisLabelFontSize: 10,
  showAxis: false,
  gridPadding: 5,
  tickLength: 5,
  tickTextPadding: 2,
  precision: 1,
  useStroke: true,
  metricPrefix: "k",
  hide: false,
};

export const defaultLineConfig = {
  trackPadding: 10,
  trackWidth: 40,
  textPadding: 5,
  textFontSize: 10,
  hide: false,
  innerRadius: 200,
  colour: "#000000",
};

export const defaultAnnotationConfig = {
  trackPadding: 5,
  trackWidth: 20,
  hide: false,
  innerRadius: 200,
  textFontSize: 10,
  textPadding: 3,
};

export const defaultGlobalConfig = {
  canvasWidth: 750,
  canvasHeight: 750,
  linkSelectedOpacity: 0.6,
  linkUnselectedOpacity: 0.1,
};
