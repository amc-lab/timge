import React, { useState } from "react";
import Box from "@mui/joy/Box";
import Slider from "@mui/joy/Slider";
import { GlobalConfig } from "@/app/types/genomes";

interface GlobalConfigProps {
  onUpdate: (newConfig: GlobalConfig) => void;
  defaultConfig: GlobalConfig;
}

export const GlobalForm: React.FC<GlobalConfigProps> = ({
  onUpdate,
  defaultConfig,
}) => {
  const [config, setConfig] = useState(defaultConfig);

  const handleConfigChange = (key: keyof GlobalConfig, value: any) => {
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
          <label>Canvas Width</label>
          <Slider
            valueLabelDisplay="auto"
            variant="solid"
            min={0}
            max={1000}
            step={10}
            value={config.canvasWidth}
            onChange={(e, value) => handleConfigChange("canvasWidth", value)}
          />
          <label>Canvas Height</label>
          <Slider
            valueLabelDisplay="auto"
            variant="solid"
            min={0}
            max={1000}
            step={10}
            value={config.canvasHeight}
            onChange={(e, value) => handleConfigChange("canvasHeight", value)}
          />
        </Box>
      </form>
    </div>
  );
};
