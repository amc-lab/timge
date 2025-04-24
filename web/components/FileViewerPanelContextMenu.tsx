import React, { useEffect, useRef, useState } from "react";
import { Menu, MenuItem } from "@mui/joy";

interface ContextMenuProps {
  visible: boolean;
  x: number;
  y: number;
  onClose: () => void;
  onCreateFolder: () => void;
  onDelete: () => void;
  onRename: () => void;
}

const ContextMenu: React.FC<ContextMenuProps> = ({
  visible,
  x,
  y,
  onClose,
  onCreateFolder,
  onDelete,
  onRename,
}) => {
  const anchorRef = useRef<HTMLDivElement>(null);
  const [anchorEl, setAnchorEl] = useState<HTMLElement | null>(null);

  useEffect(() => {
    if (visible && anchorRef.current) {
      setAnchorEl(anchorRef.current);
    }
  }, [visible, x, y]);

  useEffect(() => {
    const handleClickAway = (e: MouseEvent) => {
      if (
        anchorRef.current &&
        !anchorRef.current.contains(e.target as Node)
      ) {
        onClose();
      }
    };

    if (visible) {
      document.addEventListener("mousedown", handleClickAway);
    }

    return () => {
      document.removeEventListener("mousedown", handleClickAway);
    };
  }, [visible, onClose]);

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
            pointerEvents: "none",
          }}
        />
      )}

      <Menu
        open={visible}
        onClose={onClose}
        anchorEl={anchorEl}
        variant="soft"
        sx={{ minWidth: 180,
          fontSize: "0.8rem",
         }}
      >
        <MenuItem onClick={() => { onCreateFolder(); onClose(); }}>Create folder</MenuItem>
        <MenuItem onClick={() => { onRename(); onClose(); }}>Rename</MenuItem>
        <MenuItem onClick={() => { onDelete(); onClose(); }}>Delete</MenuItem>
      </Menu>
    </>
  );
};

export default ContextMenu;
