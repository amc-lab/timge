import React, { useState } from 'react';
import Modal from '@mui/joy/Modal';
import ModalDialog from '@mui/joy/ModalDialog';
import Typography from '@mui/joy/Typography';
import TextField from '@mui/joy/TextField';
import Button from '@mui/joy/Button';
import Stack from '@mui/joy/Stack';
import IconButton from '@mui/joy/IconButton';
import Box from '@mui/joy/Box';
import Input from '@mui/joy/Input';
import Add from '@mui/icons-material/Add';
import Remove from '@mui/icons-material/Remove';
import { FormControl, FormLabel } from '@mui/joy';
import { API_BASE_URL } from "@/app/config/env";

export default function ShapeReactivityForm({ open, onClose, onSubmit }) {
  const [conditions, setConditions] = useState([
    { name: '', files: [] }
  ]);

  const handleNameChange = (index, value) => {
    setConditions(conds => conds.map((c, i) => i === index ? { ...c, name: value } : c));
  };

  const handleFileChange = (index, event) => {
    const files = Array.from(event.target.files);
    setConditions(conds => conds.map((c, i) => i === index ? { ...c, files } : c));
  };

  const addCondition = () => {
    setConditions(conds => [...conds, { name: '', files: [] }]);
  };

  const removeCondition = (index) => {
    setConditions(conds => conds.filter((_, i) => i !== index));
  };

  const handleSubmit = () => {
    const formData = new FormData();
    conditions.forEach((cond, idx) => {
      formData.append(`condition_${idx}_name`, cond.name);
      cond.files.forEach((file, fileIdx) => {
        formData.append(`condition_${idx}_file_${fileIdx}`, file);
      });
    });
    fetch(`${API_BASE_URL}/api/timge/diff_structure/`, {
      method: 'POST',
      body: formData,
    })
    .then(response => {
      if (!response.ok) {
        throw new Error('Network response was not ok');
      }
      return response.blob();
    })
    .then(blob => {
      const url = window.URL.createObjectURL(blob);
      const a = document.createElement("a");
        a.href = url;
        a.download = "multilift.zip";
        document.body.appendChild(a);
        a.click();
        a.remove();
        window.URL.revokeObjectURL(url);
    })
    onSubmit();
  };

  return (
    <Box sx={{ 
        position: "fixed",
        top: 0,
        left: 0,
        width: "100vw",
        height: "100vh",
        backgroundColor: "rgba(0, 0, 0, 0.5)",
        zIndex: 1299, 
        display: "flex",
        justifyContent: "center",
        alignItems: "center",
        }}>
    <Modal open={open} onClose={onClose}>
      <ModalDialog
        aria-labelledby="shape-form-title"
        sx={{
          width: '60%',
        }}
      >
        <Typography id="shape-form-title" level="h4" mb={1}>
          Upload SHAPE Reactivity Data
        </Typography>
        <Typography mb={2}>
          Add conditions and upload replicate files below.
        </Typography>

        {conditions.map((cond, idx) => (
          <Box key={idx} sx={{ mb: 2, p: 2, border: '1px solid', borderColor: 'divider', borderRadius: 'md', boxShadow: 'sm' }}>
            <Stack direction="row" spacing={1} alignItems="center">
              <FormControl sx={{ flex: 1 }}>
                <FormLabel>Condition {idx + 1} Name</FormLabel>
                <Input
                  value={cond.name}
                  onChange={(e) => handleNameChange(idx, e.target.value)}
                  placeholder="Enter condition name"
                />
              </FormControl>
              <IconButton
                variant="outlined"
                color="danger"
                onClick={() => removeCondition(idx)}
                disabled={conditions.length === 1}
              >
                <Remove />
              </IconButton>
            </Stack>
            <Box component="input" type="file" multiple onChange={(e) => handleFileChange(idx, e)} sx={{ mt: 1, width: '100%' }} />
          </Box>
        ))}

        <Button
          startDecorator={<Add />}
          variant="outlined"
          onClick={addCondition}
          sx={{ mb: 2 }}
        >
          Add Condition
        </Button>

        <Stack direction="row" spacing={2} justifyContent="flex-end">
          <Button variant="plain" color="neutral" onClick={onClose}>
            Cancel
          </Button>
          <Button variant="solid" color="primary" onClick={handleSubmit}>
            Submit
          </Button>
        </Stack>
      </ModalDialog>
    </Modal>
    </Box>
  );
}
