"use client";
import React, { useEffect, useState } from "react";
import {
  Box,
  IconButton,
  Typography,
  Menu,
  MenuItem,
} from "@mui/material";
import MenuIcon from "@mui/icons-material/Menu";
import CloseIcon from "@mui/icons-material/Close";
import { useAppDispatch, useAppSelector } from "@/store/hooks";
import { deleteView, deleteConnection, deleteDependency } from "@/store/features/space/spaceSlice";

interface ParentViewProps {
  children?: React.ReactNode;
  viewConfig?: any;
  userActions?: Record<string, (...args: any[]) => void>;
  index?: number;
  ref?: React.Ref<any>;
}

const ParentView: React.FC<ParentViewProps> = ({ children, viewConfig, userActions = {}, index, ref }) => {
  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
  const open = Boolean(anchorEl);

  const dispatch = useAppDispatch();
  const space = useAppSelector((state) => state.space);

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

  const handleViewClose = () => {
    const uuid = viewConfig.uuid;
    if (uuid) {
      dispatch(deleteView(uuid));
      dispatch(deleteConnection(uuid));
      dispatch(deleteDependency(uuid));
    }
  }

  useEffect(() => {
    console.log(viewConfig);
  }, []);

  return (
    <Box
      sx={{
        display: "flex",
        justifyContent: "center",
        alignItems: "center",
        width: viewConfig.config.isMinimised ? "calc(50% - 5px)" : "calc(100% - 5px)",
        borderRadius: "3px",
        margin: "2.5px",
        flexDirection: "column",
        backgroundColor: "white",
        border: "4px solid darkblue",
      }}
      ref={ref}
    >
      <Box
        sx={{
          display: "flex",
          alignItems: "center",
          width: "100%",
          height: "2em",
          backgroundColor: "darkblue",
          px: 1,
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

        <Menu anchorEl={anchorEl} open={open} onClose={handleMenuClose}>
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

        <Box sx={{ ml: "auto" }}>
          <IconButton
            onClick={handleViewClose}
            sx={{
              color: "white",
              height: "2em",
              "&:hover": { background: "none", color: "white" },
              padding: "0px",
            }}
          >
            <CloseIcon />
          </IconButton>
        </Box>
      </Box>

      <Box
        sx={{
          backgroundColor: "white",
          borderRadius: "3px",
          width: "100%",
          minHeight: "200px",
          // padding: "5px",
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
