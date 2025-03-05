"use client";
import { useRef, useState } from "react";
import { Box, Menu, MenuItem, Button } from "@mui/joy";
import Link from "next/link";

interface DropdownMenuProps {
  label: string;
  items: { 
    text: string,
    link?: string, 
    action?: () => void,
  }[];
}

const DropdownMenu: React.FC<DropdownMenuProps> = ({ label, items }) => {

  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
  const inputRef = useRef<HTMLInputElement>(null);
  const handleMouseEnter = (event: React.MouseEvent<HTMLDivElement>) => {
    setAnchorEl(event.currentTarget);
  };

  const handleMouseLeave = () => {
    setAnchorEl(null);
  };

  const triggerFileSelect = () => {
    if (inputRef.current) {
      inputRef.current.click();
    }
  };

  return (
    <Box
      onMouseEnter={handleMouseEnter}
      onMouseLeave={handleMouseLeave}
      sx={{ position: "relative", height: "100%" }} // Ensure full height
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
          minWidth: "4em"
        }}
      >
        {label}
      </Button>
      <Menu
        anchorEl={anchorEl}
        open={Boolean(anchorEl)}
        onClose={handleMouseLeave}
        placement="bottom-start"
      >
        {items.map((item, index) => (
          <MenuItem key={index} onClick={handleMouseLeave}>
            {item.action ? (
                <Box
                  onClick={item.action}
                  sx={{
                    background: "none",
                    color: "black",
                  }}
                >
                  {item.text}
                </Box>
              ) : (
              <Link href={item.link}>{item.text}</Link>
              )
            }
          </MenuItem>
        ))}
      </Menu>
    </Box>
  );
};

export default DropdownMenu;
