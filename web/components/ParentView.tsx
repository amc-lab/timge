"use client";
import React from "react";
import { Box, IconButton, Typography } from "@mui/material";
import MenuIcon from "@mui/icons-material/Menu";

interface ParentViewProps {
  children?: React.ReactNode;
  viewConfig?: any;
}

const ParentView: React.FC<ParentViewProps> = ({ children, viewConfig }) => {
    console.log("ParentView props", viewConfig);
  return (
    <Box
      sx={{
        display: "flex",
        justifyContent: "center",
        alignItems: "center",
        width: "calc(100% - 5px)",
        borderRadius: "3px",
        margin: "2.5px",
        flexDirection: "column",
        backgroundColor: "white",
        border: "4px solid darkblue",
      }}
    >
      <Box
        sx={{
          display: "flex",
          alignItems: "center",
          width: "100%",
          height: "2em",
          backgroundColor: "darkblue",
        }}
      >
        <IconButton
          sx={{
            color: "white",
            height: "2em",
            "&:hover": { background: "none", color: "white" },
          }}
        >
          <MenuIcon />
        </IconButton>
        <Typography
            sx={{
                color: "white",
                fontSize: "0.9em",
                fontWeight: "bold",
                marginLeft: "10px",
                }}
        >
            <strong>{viewConfig.title}</strong>
        </Typography>
      </Box>

      <Box
        sx={{
          backgroundColor: "white",
          borderRadius: "3px",
          width: "650px",
          minHeight: "200px",
          padding: "5px",
          display: "flex",
          justifyContent: "center",
          alignItems: "center",
        }}
      >
        {children}
      </Box>
    </Box>
  );
};

export default ParentView;
