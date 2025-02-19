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

export const defaultGlobalConfig = {
  canvasWidth: 750,
  canvasHeight: 750,
};
