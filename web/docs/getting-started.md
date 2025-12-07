# Getting Started

Use this quick-start to set expectations before diving into the specialised
views. TIMGE stitches together track management, multiple visualisations, and
cross-view interactions; understanding the core loop up front makes feature
pages like Multilift or vRISE much easier to follow.

## 1. Prerequisites

- **Node.js / npm** – the `web/` workspace is a Next.js app. Use Node 18+ and
  run `npm install` inside `web/`.
- **Django API** – backend services (tile generation, Multilift, file uploads)
  come from the `backend/` folder. Follow its README to bring up the server,
  then point `NEXT_PUBLIC_DJANGO_HOST` to it.
- **Data directory** – TIMGE expects uploaded tracks to live under
  `<TRACK_ROOT_DIR>/<uuid>/`. Create a demo UUID via the UI or seed the folder
  manually with BED, BEDPE, BEDGraph, and FASTA/Fasta index pairs.

## 2. Launching the Workspace

```bash
cd web
npm run dev
```

Open `http://localhost:3000` and upload a few tracks from the sidebar. The
uploaded filenames populate `space.dataFiles`, which all views reference when
asking for segments or contact maps.

## 3. Core Concepts

- **Views** – linear, circos, and map views can be added via the header
  dropdown. Each view maintains its own config but can subscribe to dependency
  updates from other views.
- **Track selector** – uniform modal for picking track files; supports clicking
  on list items rather than checkboxes for faster selection.
- **Cross-view events** – propagated through a lightweight pub/sub helper
  (`crossViewEvents.ts`). This powers actions like “Generate heatmap” from the
  circos context menu.

## 4. Recommended Onboarding Flow

1. Upload (or reuse) a small FASTA + BED/BEDGraph set.
2. Add a linear view to inspect tracks in IGV.js.
3. Add a circos view using the same inputs and experiment with selections.
4. Add a map view, choose segments, and render tiles.
5. Try Multilift to perform coordinate conversions and see how its output feeds
   back into other views.

## 5. One-command Local Setup

For contributors who just want the stack running locally, use the helper script
at the repo root:

```bash
./run-local.sh
```

The script performs the following:

1. Installs backend deps via Poetry and runs migrations.
2. Installs frontend deps via npm, builds the production bundle, and starts the
   Next.js server.
3. Exposes sensible defaults (`NEXT_PUBLIC_DJANGO_HOST=http://127.0.0.1:8000`,
   `NEXT_PUBLIC_FILE_HOST=http://127.0.0.1:8000/uploads/`, uploads stored in
   `<repo>/uploads`).
4. Creates a symlink from `web/public/uploads` to the upload directory so IGV
   and other viewers can resolve file URLs.
5. Traps `Ctrl+C` to stop both servers cleanly.

Override ports or directories if needed:

```bash
FRONTEND_PORT=4000 BACKEND_PORT=9000 TRACK_ROOT_DIR=/tmp/timge_uploads ./run-local.sh
```

Make sure `poetry` and `npm` are installed before running the script.

## 5. Troubleshooting

- **404s on tile endpoints** – confirm the Django host env var and ensure the
  server has read access to the track directory.
- **Hydration mismatches** – most UI state hydrates from Redux + localStorage.
  Clear localStorage (`localStorage.removeItem("timge_space_state")`) if you see
  inconsistent views on refresh.
- **Large uploads** – the UI does not currently chunk uploads; keep demo files
  small or extend the API to support streaming uploads.

When you are comfortable with this loop, continue with dedicated guides:
[Multilift](./multilift.md) and [vRISE Map View](./vrise.md).
