#!/usr/bin/env bash

set -Eeuo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FRONTEND_DIR="$ROOT_DIR/web"
BACKEND_DIR="$ROOT_DIR/backend"

FRONTEND_PORT="${FRONTEND_PORT:-3000}"
BACKEND_PORT="${BACKEND_PORT:-8000}"
TRACK_ROOT_DIR="${TRACK_ROOT_DIR:-$ROOT_DIR/uploads}"
NEXT_PUBLIC_DJANGO_HOST="${NEXT_PUBLIC_DJANGO_HOST:-http://127.0.0.1:${BACKEND_PORT}}"
NEXT_PUBLIC_FILE_HOST="${NEXT_PUBLIC_FILE_HOST:-${NEXT_PUBLIC_DJANGO_HOST}/uploads/}"

command -v npm >/dev/null 2>&1 || {
  echo "npm is required but not installed. Please install Node.js 18+." >&2
  exit 1
}

command -v poetry >/dev/null 2>&1 || {
  echo "Poetry is required but not installed. See https://python-poetry.org/docs/." >&2
  exit 1
}

mkdir -p "$TRACK_ROOT_DIR"

cleanup() {
  echo ""
  echo "Stopping local TIMGE..."
  if [[ -n "${FRONTEND_PID:-}" ]] && ps -p "$FRONTEND_PID" >/dev/null 2>&1; then
    kill "$FRONTEND_PID" >/dev/null 2>&1 || true
    wait "$FRONTEND_PID" 2>/dev/null || true
  fi
  if [[ -n "${BACKEND_PID:-}" ]] && ps -p "$BACKEND_PID" >/dev/null 2>&1; then
    kill "$BACKEND_PID" >/dev/null 2>&1 || true
    wait "$BACKEND_PID" 2>/dev/null || true
  fi
  echo "Done."
}
trap cleanup EXIT INT TERM

echo "Installing backend dependencies with Poetry..."
pushd "$BACKEND_DIR" >/dev/null
poetry install
poetry run python manage.py migrate
echo "Starting Django backend on port ${BACKEND_PORT}..."
TRACK_ROOT_DIR="$TRACK_ROOT_DIR" poetry run python manage.py runserver "0.0.0.0:${BACKEND_PORT}" &
BACKEND_PID=$!
popd >/dev/null

PUBLIC_UPLOADS_DIR="$FRONTEND_DIR/public/uploads"
if [[ ! -e "$PUBLIC_UPLOADS_DIR" ]]; then
  mkdir -p "$(dirname "$PUBLIC_UPLOADS_DIR")"
  ln -s "$TRACK_ROOT_DIR" "$PUBLIC_UPLOADS_DIR"
fi

echo "Installing frontend dependencies..."
pushd "$FRONTEND_DIR" >/dev/null
npm install

echo "Building production bundle..."
NEXT_PUBLIC_DJANGO_HOST="$NEXT_PUBLIC_DJANGO_HOST" \
NEXT_PUBLIC_FILE_HOST="$NEXT_PUBLIC_FILE_HOST" \
npm run build

echo "Starting Next.js server on port ${FRONTEND_PORT}..."
NEXT_PUBLIC_DJANGO_HOST="$NEXT_PUBLIC_DJANGO_HOST" \
NEXT_PUBLIC_FILE_HOST="$NEXT_PUBLIC_FILE_HOST" \
npm run start -- --hostname 0.0.0.0 --port "$FRONTEND_PORT" &
FRONTEND_PID=$!
popd >/dev/null

echo ""
echo "TIMGE is live:"
echo "  Frontend:  http://localhost:${FRONTEND_PORT}"
echo "  Backend:   ${NEXT_PUBLIC_DJANGO_HOST}"
echo "  Uploads:   ${TRACK_ROOT_DIR}"
echo ""
echo "Press Ctrl+C to stop both services."

wait "$FRONTEND_PID"
