'use client';

import { Input, Typography, Box } from '@mui/joy';

export default function EditableLabel({
  value,
  isEditing,
  onChange,
  onSave,
}: {
  value: string;
  isEditing: boolean;
  onChange: (val: string) => void;
  onSave: () => void;
}) {
  return (
    <Box
        sx={{
            display: "flex",
            alignItems: "center",
            height: "100%",
        }}
    >
      {isEditing ? (
        <Input
          value={value}
          onChange={(e) => onChange(e.target.value)}
          onKeyDown={(e) => {
            if (e.key === 'Enter') {
              e.preventDefault();
              onSave();
            }
          }}
          autoFocus
          variant="plain"
          size="sm"
          sx={{
            border: "none",
            borderRadius: "3px",
            '--Input-minHeight': '1.75em',
            '--Input-paddingInline': '0.25em',
            '--Input-paddingBlock': '0',
            '&:focus-within': {
                outline: 'none',
                boxShadow: 'none',
              },
              '& input': {
                outline: 'none',
                boxShadow: 'none',
              },
          }}
        />
      ) : (
        <Typography
          level="body-sm"
          sx={{
            color: "white",
            fontSize: "0.9em",
            fontWeight: "bold",
            marginLeft: "10px",
            lineHeight: "2em",
          }}
        >
          {value}
        </Typography>
      )}
    </Box>
  );
}
