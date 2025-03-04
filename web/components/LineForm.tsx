import React from "react";
import Box from "@mui/joy/Box";
import Slider from "@mui/joy/Slider";
import { LineConfig } from "@/app/types/genomes";

interface LineFormProps {
  onUpdate: (newConfig: LineConfig) => void;
  defaultConfig: LineConfig;
}

export const LineForm: React.FC<LineFormProps> = ({
  onUpdate,
  defaultConfig,
}) => {
  const [config, setConfig] = React.useState(defaultConfig);

  const handleConfigChange = (key: keyof LineConfig, value: any) => {
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
          <label>Track Width</label>
          <Slider
            valueLabelDisplay="auto"
            variant="solid"
            min={20}
            max={50}
            step={2}
            value={defaultConfig.trackWidth}
            onChange={(e, value) => handleConfigChange("trackWidth", value)}
          />
          <label>Track Padding</label>
          <Slider
            valueLabelDisplay="auto"
            variant="solid"
            min={0}
            max={20}
            step={2}
            value={defaultConfig.trackPadding}
            onChange={(e, value) => handleConfigChange("trackPadding", value)}
          />
          <label>Colour</label>
            <input
                type="color"
                value={defaultConfig.colour}
                onChange={(e) => handleConfigChange("colour", e.target.value)}
            />
        </Box>
      </form>
    </div>
  );
};
