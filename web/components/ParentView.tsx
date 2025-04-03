"use client";
import React, { useState } from "react";
import {
  Box,
  IconButton,
  Typography,
  Menu,
  MenuItem,
} from "@mui/material";
import MenuIcon from "@mui/icons-material/Menu";

interface ParentViewProps {
  children?: React.ReactNode;
  viewConfig?: any;
  userActions?: Record<string, (...args: any[]) => void>;
}

const ParentView: React.FC<ParentViewProps> = ({ children, viewConfig, userActions = {} }) => {
  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
  const open = Boolean(anchorEl);

  const handleMenuClick = (event: React.MouseEvent<HTMLElement>) => {
    setAnchorEl(event.currentTarget);
  };

  const handleMenuClose = () => {
    setAnchorEl(null);
  };

  const handleActionClick = (actionKey: string) => {
    userActions[actionKey]?.();
    handleMenuClose();
  };

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
          onClick={handleMenuClick}
          sx={{
            color: "white",
            height: "2em",
            "&:hover": { background: "none", color: "white" },
          }}
        >
          <MenuIcon />
        </IconButton>
        <Menu
          anchorEl={anchorEl}
          open={open}
          onClose={handleMenuClose}
        >
          {Object.keys(userActions).map((key) => (
            <MenuItem key={key} onClick={() => handleActionClick(key)}>
              {key}
            </MenuItem>
          ))}
        </Menu>
        <Typography
          sx={{
            color: "white",
            fontSize: "0.9em",
            fontWeight: "bold",
            marginLeft: "10px",
          }}
        >
          <strong>{viewConfig?.title}</strong>
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
