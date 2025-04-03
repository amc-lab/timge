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

const LoginPage: React.FC = () => {
  const [email, setEmail] = useState("");
  const [password, setPassword] = useState("");
  const [error, setError] = useState("");

  const handleLogin = (e: React.FormEvent) => {
    e.preventDefault();

    // Example validation or API call
    if (!email || !password) {
      setError("Please fill in all fields.");
      return;
    }

    setError("");
    console.log({ email, password });
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
          Login
        </Typography>
        <Divider sx={{ mb: 2 }} />

        {error && (
          <Alert color="danger" variant="soft" sx={{ mb: 2 }}>
            {error}
          </Alert>
        )}

        <form onSubmit={handleLogin}>
          <Stack spacing={2}>
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

            <Typography level="body-sm" textAlign="center">
              New here?{" "}
              <Link href="/register" passHref legacyBehavior>
                <Typography
                  component="a"
                  level="body-sm"
                  sx={{ textDecoration: "underline", cursor: "pointer" }}
                  color="primary"
                >
                  Register
                </Typography>
              </Link>
            </Typography>

            <Button type="submit" color="primary">
              Login
            </Button>
          </Stack>
        </form>
      </Card>
    </Box>
  );
};

export default LoginPage;
