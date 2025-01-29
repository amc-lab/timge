import React from "react";
import { BarConfig } from "@/app/types/genomes";
import Box from "@mui/joy/Box";
import Slider from "@mui/joy/Slider";

interface BarFormProps {
  onUpdate: (newConfig: BarConfig) => void;
  defaultConfig: BarConfig;
}

export const BarForm: React.FC<BarFormProps> = ({
  onUpdate,
  defaultConfig,
}) => {
  const [config, setConfig] = React.useState(defaultConfig);

  const handleConfigChange = (key: keyof BarConfig, value: any) => {
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
            min={0}
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
        </Box>
      </form>
    </div>
  );
};
