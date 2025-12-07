"use client";
import { useRef, useState } from "react";
import { Box, Menu, MenuItem, Button } from "@mui/joy";
import Link from "next/link";

interface DropdownMenuProps {
  label: string;
  items: {
    text: string;
    link?: string;
    action?: () => void;
  }[];
}

const DropdownMenu: React.FC<DropdownMenuProps> = ({ label, items }) => {
  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
  const [menuOpen, setMenuOpen] = useState(false);
  const closeTimeoutRef = useRef<NodeJS.Timeout | null>(null);

  const handleMouseEnter = (event: React.MouseEvent<HTMLDivElement>) => {
    if (closeTimeoutRef.current) {
      clearTimeout(closeTimeoutRef.current);
    }
    setAnchorEl(event.currentTarget);
    setMenuOpen(true);
  };

  const handleMouseLeave = () => {
    closeTimeoutRef.current = setTimeout(() => {
      setMenuOpen(false);
      setAnchorEl(null);
    }, 100); // Delay in ms
  };

  return (
    <Box
      onMouseEnter={handleMouseEnter}
      onMouseLeave={handleMouseLeave}
      sx={{ position: "relative", height: "100%" }}
    >
      <Button
        variant="solid"
        sx={{
          backgroundColor: "black",
          color: "white",
          height: "100%",
          display: "flex",
          alignItems: "center",
          px: 3,
          "&:hover": { backgroundColor: "#333" },
          minWidth: "4em",
        }}
      >
        {label}
      </Button>
      <Menu
        anchorEl={anchorEl}
        open={menuOpen}
        onClose={() => setMenuOpen(false)}
        placement="bottom-start"
        onMouseEnter={() => {
          if (closeTimeoutRef.current) clearTimeout(closeTimeoutRef.current);
        }}
        onMouseLeave={handleMouseLeave}
      >
        {items.map((item, index) => (
          <MenuItem key={index} 
            onClick={() => {
              setMenuOpen(false);
              if (item.action) item.action();
              }}>
            {item.action ? (
              <Box sx={{ background: "none", color: "black" }}>
                {item.text}
              </Box>
            ) : (
              <Link href={item.link || "#"}>{item.text}</Link>
            )}
          </MenuItem>
        ))}
      </Menu>
    </Box>
  );
};

export default DropdownMenu;
