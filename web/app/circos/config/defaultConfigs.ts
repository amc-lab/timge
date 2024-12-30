export const defaultAssemblyConfig = {
    segmentPadding: 0.02,
    axisLabelFontSize: 10,
    showAxis: true,
    segmentInnerRadius: 200,
    segmentOuterRadius: 250,
    segmentGridPadding: 5,
    tickLength: 5,
    tickTextPadding: 2,
    precision: 1,
    useStroke: true,
    metricPrefix: "M"
};

export const defaultChordConfig = {
    chordPadding: 10,
    opacity: 0.3,
    color: "blue",
    useStroke: true,
    outerRadius: defaultAssemblyConfig.segmentInnerRadius
};

export const defaultGlobalConfig = {
    canvasWidth: 650,
    canvasHeight: 650,
};