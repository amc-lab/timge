import React, { useState } from "react";
import Checkbox from "@mui/joy/Checkbox";
import Box from "@mui/joy/Box";
import Slider from "@mui/joy/Slider";
import { RingConfig } from "@/app/types/genomes";
import Select from "@mui/joy/Select";
import Option from "@mui/joy/Option";

interface RingFormProps {
  onUpdate: (newConfig: RingConfig) => void;
  defaultConfig: RingConfig;
}

export const RingForm: React.FC<RingFormProps> = ({
  onUpdate,
  defaultConfig,
}) => {
  const [config, setConfig] = useState(defaultConfig);

  const handleConfigChange = (key: keyof RingConfig, value: any) => {
    const updatedConfig = { ...config, [key]: value };
    setConfig(updatedConfig);
    onUpdate(updatedConfig);
  };

  return (
    <div>
      <form>
        <Box
          sx={{
            display: "grid",
            gridTemplateColumns: "20% 80%",
            gap: 1,
            alignItems: "center",
            width: "100%",
            padding: 1,
            fontSize: "0.65rem",
          }}
        >
          {/* <label>Segment Padding</label>
          <Slider
            valueLabelDisplay="auto"
            variant="solid"
            min={0}
            max={0.1}
            step={0.01}
            value={config.segmentPadding}
            onChange={(e, value) => handleConfigChange("segmentPadding", value)}
          /> */}
          <label>Axis Font Size</label>
          <Slider
            valueLabelDisplay="auto"
            variant="solid"
            min={6}
            max={16}
            step={1}
            value={config.axisLabelFontSize}
            onChange={(e, value) =>
              handleConfigChange("axisLabelFontSize", value)
            }
          />
          <label>Show Axis</label>
          <Checkbox
            label=""
            checked={config.showAxis}
            onChange={(e) => handleConfigChange("showAxis", e.target.checked)}
          />
          <label>Track Width</label>
          <Slider
            valueLabelDisplay="auto"
            variant="solid"
            min={20}
            max={50}
            step={2}
            value={config.trackWidth}
            onChange={(e, value) =>
              handleConfigChange("trackWidth", value)
            }
          />
          <label>Track Padding</label>
          <Slider
            valueLabelDisplay="auto"
            variant="solid"
            min={5}
            max={25}
            step={1}
            value={config.trackPadding}
            onChange={(e, value) =>
              handleConfigChange("trackPadding", value)
            }
          />
          <label>Grid Padding</label>
          <Slider
            valueLabelDisplay="auto"
            variant="solid"
            min={0}
            max={5}
            step={1}
            value={config.gridPadding}
            onChange={(e, value) =>
              handleConfigChange("gridPadding", value)
            }
          />
          <label>Tick Length</label>
          <Slider
            valueLabelDisplay="auto"
            variant="solid"
            min={0}
            max={10}
            step={1}
            value={config.tickLength}
            onChange={(e, value) => handleConfigChange("tickLength", value)}
          />
          <label>Text Padding</label>
          <Slider
            valueLabelDisplay="auto"
            variant="solid"
            min={0}
            max={10}
            step={1}
            value={config.tickTextPadding}
            onChange={(e, value) =>
              handleConfigChange("tickTextPadding", value)
            }
          />
          <label>Precision</label>
          <Slider
            valueLabelDisplay="auto"
            variant="solid"
            min={0}
            max={2}
            step={1}
            value={config.precision}
            onChange={(e, value) => handleConfigChange("precision", value)}
          />
          <label>Use Stroke</label>
          <Checkbox
            label=""
            checked={config.useStroke}
            onChange={(e) => handleConfigChange("useStroke", e.target.checked)}
          />
          <label>Metric Prefix</label>
          <Select
            defaultValue="M"
            onChange={(e, value) => handleConfigChange("metricPrefix", value)}
            size="sm"
          >
            <Option value="k">k</Option>
            <Option value="M">M</Option>
            <Option value="G">G</Option>
            <Option value="None">None</Option>
          </Select>
          <label>Hide</label>
          <Checkbox
            label=""
            checked={config.hide}
            onChange={(e) => handleConfigChange("hide", e.target.checked)}
          />
        </Box>
      </form>
    </div>
  );
};
