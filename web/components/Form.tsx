import React, { useState } from "react";

export default function AssemblyForm({ onUpdate }: { 
    onUpdate: (
        segmentPadding: number,
        axisLabelFontSize: number, 
        showAxis: boolean,
        segmentInnerRadius: number,
        segmentOuterRadius: number,
        segmentGridPadding: number,
        canvasWidth: number,
        canvasHeight: number,
        tickLength: number,
        tickTextPadding: number,
        precision: number
    ) => void 
}) {
  const [segmentPadding, setSegmentPadding] = useState(0.02);
  const [axisLabelFontSize, setAxisLabelFontSize] = useState(10);
  const [showAxis, setShowAxis] = useState(true);
  const [segmentInnerRadius, setSegmentInnerRadius] = useState(200);
  const [segmentOuterRadius, setSegmentOuterRadius] = useState(250);
  const [segmentGridPadding, setSegmentGridPadding] = useState(3);
  const [canvasWidth, setCanvasWidth] = useState(600);
  const [canvasHeight, setCanvasHeight] = useState(600);
  const [tickLength, setTickLength] = useState(3);
  const [tickTextPadding, setTickTextPadding] = useState(3);
  const [precision, setPrecision] = useState(1);

    const onFormUpdate = () => {
        if (segmentInnerRadius >= segmentOuterRadius) {
            alert("Segment inner radius must be less than segment outer radius");
            setSegmentInnerRadius(200);
            setSegmentOuterRadius(250);
            return;
        }
        onUpdate(segmentPadding, axisLabelFontSize, showAxis, segmentInnerRadius, segmentOuterRadius, segmentGridPadding, canvasWidth, canvasHeight, tickLength, tickTextPadding, precision);
    };

  return (
    <div>
        <form>
            <label>
                Segment Padding</label>
                <input type="range" min="0" max="0.1" step="0.01" value={segmentPadding} onChange={(e) => setSegmentPadding(parseFloat(e.target.value))} />
            <br />
            <label>
                Axis Label Font Size</label>
                <input type="range" min="6" max="16" step="1" value={axisLabelFontSize} onChange={(e) => setAxisLabelFontSize(parseFloat(e.target.value))} />
            <br />
            <label>
                Show Axis
                <input type="checkbox" checked={showAxis} onChange={(e) => setShowAxis(e.target.checked)} />
            </label>
            <br />
            <label>
                Segment Inner Radius</label>
                <input type="range" min="0" max="500" step="10" value={segmentInnerRadius} onChange={(e) => setSegmentInnerRadius(parseFloat(e.target.value))} />
            <br />
            <label>
                Segment Outer Radius</label>
                <input type="range" min="0" max="500" step="10" value={segmentOuterRadius} onChange={(e) => setSegmentOuterRadius(parseFloat(e.target.value))} />
            <br />
            <label>
                Segment Grid Padding</label>
                <input type="range" min="1" max="5" step="1" value={segmentGridPadding} onChange={(e) => setSegmentGridPadding(parseFloat(e.target.value))} />
            <br />
            <br />
            <label>
                Canvas Width</label>
                <input type="range" min="0" max="1000" step="10" value={canvasWidth} onChange={(e) => setCanvasWidth(parseFloat(e.target.value))} />
            <br />
            <label>
                Canvas Height</label>
                <input type="range" min="0" max="1000" step="10" value={canvasHeight} onChange={(e) => setCanvasHeight(parseFloat(e.target.value))} />
            <br />
            <label>
                Tick Length</label>
                <input type="range" min="0" max="10" step="1" value={tickLength} onChange={(e) => setTickLength(parseFloat(e.target.value))} />
            <br />
            <label>
                Tick Text Padding</label>
                <input type="range" min="0" max="10" step="1" value={tickTextPadding} onChange={(e) => setTickTextPadding(parseFloat(e.target.value))} />
            <br />
            <label>
                Precision</label>
                <input type="range" min="0" max="2" step="1" value={precision} onChange={(e) => setPrecision(parseFloat(e.target.value))} />
            <br />
        </form>
      <button onClick={onFormUpdate}>Update</button>
    </div>
  );
}