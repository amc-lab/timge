import React, { useState } from "react";
import { AssemblyConfig } from "@/app/types/genomes";

interface AssemblyFormProps {
    onUpdate: (newConfig: AssemblyConfig) => void;
    defaultConfig: AssemblyConfig;
}

export const AssemblyForm: React.FC<AssemblyFormProps> = ({ onUpdate, defaultConfig }) => {
    const [config, setConfig] = useState(defaultConfig);

    const onFormUpdate = () => {
        if (config.segmentInnerRadius >= config.segmentOuterRadius) {
            alert("Segment inner radius must be less than segment outer radius");
            return;
        }
        onUpdate(config);
    };

  return (
    <div>
        <form>
            <label>
                Segment Padding</label>
                <input type="range" min="0" max="0.1" step="0.01" value={config.segmentPadding} onChange={(e) => setConfig({...config, segmentPadding: parseFloat(e.target.value)})} />
            <br />
            <label>
                Axis Label Font Size</label>
                <input type="range" min="6" max="16" step="1" value={config.axisLabelFontSize} onChange={(e) => setConfig({...config, axisLabelFontSize: parseFloat(e.target.value)})} />
            <br />
            <label>
                Show Axis
                <input type="checkbox" checked={config.showAxis} onChange={(e) => setConfig({...config, showAxis: e.target.checked})} />
            </label>
            <br />
            <label>
                Segment Inner Radius</label>
                <input type="range" min="0" max="500" step="10" value={config.segmentInnerRadius} onChange={(e) => setConfig({...config, segmentInnerRadius: parseFloat(e.target.value)})} />
            <br />
            <label>
                Segment Outer Radius</label>
                <input type="range" min="0" max="500" step="10" value={config.segmentOuterRadius} onChange={(e) => setConfig({...config, segmentOuterRadius: parseFloat(e.target.value)})} />
            <label>
                Segment Grid Padding</label>
                <input type="range" min="1" max="5" step="1" value={config.segmentGridPadding} onChange={(e) => setConfig({...config, segmentGridPadding: parseFloat(e.target.value)})} />
            <br />
            <br />
            <label>
                Canvas Width</label>
                <input type="range" min="0" max="1000" step="10" value={config.canvasWidth} onChange={(e) => setConfig({...config, canvasWidth: parseFloat(e.target.value)})} />
            <label>
                Canvas Height</label>
                <input type="range" min="0" max="1000" step="10" value={config.canvasHeight} onChange={(e) => setConfig({...config, canvasHeight: parseFloat(e.target.value)})} />
            <br />
            <label>
                Tick Length</label>
                <input type="range" min="0" max="10" step="1" value={config.tickLength} onChange={(e) => setConfig({...config, tickLength: parseFloat(e.target.value)})} />
            <br />
            <label>
                Tick Text Padding</label>
                <input type="range" min="0" max="10" step="1" value={config.tickTextPadding} onChange={(e) => setConfig({...config, tickTextPadding: parseFloat(e.target.value)})} />
            <br />
            <label>
                Precision</label>
                <input type="range" min="0" max="2" step="1" value={config.precision} onChange={(e) => setConfig({...config, precision: parseFloat(e.target.value)})} />
            <br />
            <label>
                Use Stroke
                <input type="checkbox" checked={config.useStroke} onChange={(e) => setConfig({...config, useStroke: e.target.checked})} />
            </label>
        </form>
      <button onClick={onFormUpdate}>Update</button>
    </div>
  );
}