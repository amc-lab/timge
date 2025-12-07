import React from "react";
import { Box, Button, Divider, Sheet, Typography } from "@mui/joy";

interface TrackSelectorModalProps {
  title: string;
  description?: string;
  onClose: () => void;
  onConfirm: () => void;
  confirmLabel?: string;
  confirmDisabled?: boolean;
  children: React.ReactNode;
}

const TrackSelectorModal = ({
  title,
  description,
  onClose,
  onConfirm,
  confirmLabel = "Confirm Selection",
  confirmDisabled = false,
  children,
}: TrackSelectorModalProps) => {
  return (
    <Box
      sx={{
        position: "fixed",
        top: 0,
        left: 0,
        width: "100vw",
        height: "100vh",
        backgroundColor: "rgba(0, 0, 0, 0.6)",
        display: "flex",
        justifyContent: "center",
        alignItems: "center",
        zIndex: 1300,
      }}
    >
      <Sheet
        sx={{
          width: "70vw",
          maxWidth: "1000px",
          maxHeight: "80vh",
          backgroundColor: "white",
          borderRadius: "10px",
          padding: 2.5,
          display: "flex",
          flexDirection: "column",
          gap: 1,
          boxShadow: "0 18px 45px rgba(15, 23, 42, 0.35)",
        }}
      >
        <Box sx={{ display: "flex", flexDirection: "column", gap: 0.5 }}>
          <Typography level="title-lg" sx={{ lineHeight: 1.2 }}>
            {title}
          </Typography>
          {description && (
            <Typography level="body-sm" color="neutral">
              {description}
            </Typography>
          )}
        </Box>
        <Divider />
        <Box sx={{ flexGrow: 1, overflowY: "auto", pr: 0.5 }}>{children}</Box>
        <Divider />
        <Box display="flex" justifyContent="flex-end" gap={1.5}>
          <Button variant="outlined" color="danger" onClick={onClose}>
            Cancel
          </Button>
          <Button
            variant="solid"
            onClick={onConfirm}
            disabled={confirmDisabled}
          >
            {confirmLabel}
          </Button>
        </Box>
      </Sheet>
    </Box>
  );
};

export default TrackSelectorModal;
