import React, { useEffect, useRef, useState } from "react";
import { Box, Menu, MenuItem } from "@mui/joy";

export interface ContextMenuOption {
  label: string;
  onClick: () => void;
  disabled?: boolean;
  icon?: React.ReactNode;
}

interface ContextMenuProps {
  visible: boolean;
  x: number;
  y: number;
  onClose: () => void;
  options: ContextMenuOption[];
}

import { ClickAwayListener } from "@mui/base";

const ContextMenu: React.FC<ContextMenuProps> = ({
  visible,
  x,
  y,
  onClose,
  options,
}) => {
  const anchorRef = useRef<HTMLDivElement>(null);
  const [anchorEl, setAnchorEl] = useState<HTMLElement | null>(null);

  useEffect(() => {
    if (visible && anchorRef.current) {
      setAnchorEl(anchorRef.current);
    }
  }, [visible, x, y]);

  return (
    <>
      {visible && (
        <div
          ref={anchorRef}
          style={{
            position: "fixed",
            top: y,
            left: x,
            width: 0,
            height: 0,
          }}
        />
      )}

      <ClickAwayListener onClickAway={onClose}>
        <Menu
          open={visible}
          onClose={onClose}
          anchorEl={anchorEl}
          placement="right-start"
          disablePortal
          variant="soft"
          sx={{ minWidth: 180, fontSize: "0.8rem" }}
        >
          {options.map((opt, i) => (
            <MenuItem
              key={i}
              onClick={() => {
                opt.onClick();
                onClose();
              }}
              disabled={opt.disabled}
              sx={{ display: "flex", alignItems: "center", gap: ".5em" }}
            >
              {opt.icon}
              {opt.label}
            </MenuItem>
          ))}
        </Menu>
      </ClickAwayListener>
    </>
  );
};


export default ContextMenu;
