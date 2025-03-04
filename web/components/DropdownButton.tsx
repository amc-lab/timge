"use client";
import { Box, Button } from "@mui/joy";
import Link from "next/link";

interface DropdownButtonProps {
  label: string;
  link: string;   
}

const DropdownButton: React.FC<DropdownButtonProps> = ({ label, link}) => {
  return (
    <Box
      sx={{ position: "relative", height: "100%" }} // Ensure full height
    >
    <Link href={link}>
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
    </Link>
    </Box>
  );
};

export default DropdownButton;
