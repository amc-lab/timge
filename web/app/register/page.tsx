"use client";

import React, { useState } from "react";
import Link from "next/link";
import {
  Box,
  Button,
  FormControl,
  FormLabel,
  Input,
  Typography,
  Stack,
  Alert,
  Divider,
  Card,
} from "@mui/joy";

const RegisterPage: React.FC = () => {
  const [username, setUsername] = useState("");
  const [email, setEmail] = useState("");
  const [password, setPassword] = useState("");
  const [confirmPassword, setConfirmPassword] = useState("");
  const [error, setError] = useState("");

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();

    if (password !== confirmPassword) {
      setError("Passwords do not match.");
      return;
    }

    // Reset error and proceed (replace with API call if needed)
    setError("");
    console.log({ username, email, password });
  };

  return (
    <Box
      display="flex"
      justifyContent="center"
      alignItems="center"
      minHeight="100vh"
      sx={{ backgroundColor: "#f8f9fa" }}
    >
      <Card
        variant="outlined"
        sx={{
          padding: 4,
          minWidth: 400,
          borderRadius: "md",
          boxShadow: "lg",
        }}
      >
        <Typography level="h4" textAlign="center" mb={2}>
          Registration
        </Typography>
        <Divider sx={{ mb: 2 }} />

        {error && (
          <Alert color="danger" variant="soft" sx={{ mb: 2 }}>
            {error}
          </Alert>
        )}

        <form onSubmit={handleSubmit}>
          <Stack spacing={2}>
            <FormControl>
              <FormLabel>Username</FormLabel>
              <Input
                value={username}
                onChange={(e) => setUsername(e.target.value)}
                required
              />
            </FormControl>

            <FormControl>
              <FormLabel>Email</FormLabel>
              <Input
                type="email"
                value={email}
                onChange={(e) => setEmail(e.target.value)}
                required
              />
            </FormControl>

            <FormControl>
              <FormLabel>Password</FormLabel>
              <Input
                type="password"
                value={password}
                onChange={(e) => setPassword(e.target.value)}
                required
              />
            </FormControl>

            <FormControl>
              <FormLabel>Confirm Password</FormLabel>
              <Input
                type="password"
                value={confirmPassword}
                onChange={(e) => setConfirmPassword(e.target.value)}
                required
              />
            </FormControl>

            <Typography level="body-sm" textAlign="center">
              Already have an account?{" "}
              <Link href="/login" passHref legacyBehavior>
                <Typography
                  component="a"
                  level="body-sm"
                  sx={{ textDecoration: "underline", cursor: "pointer" }}
                  color="primary"
                >
                  Login
                </Typography>
              </Link>
            </Typography>

            <Button type="submit" color="primary">
              Create Account
            </Button>
          </Stack>
        </form>
      </Card>
    </Box>
  );
};

export default RegisterPage;
