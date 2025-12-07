# vRISE (Interactive Map View)

vRISE provides the tile-based contact map renderer inside TIMGE. It focuses on
rendering Hi-C style matrices efficiently, caching tiles, and coordinating with
other views via the event bus.

## Purpose

- Explore interaction matrices at multiple resolutions without re-rendering the
  entire dataset.
- Compare loci selected in the circos or linear view by piping the selections
  into the map view.
- Export publication-ready PNGs from the canvas/SVG hybrid renderer.

## Workflow

1. **Select tracks** – open the track picker and choose the reference FASTA and
   contact track (BEDPE/BAM-derived matrix). These choices persist via Redux +
   localStorage so refreshes keep working data.
2. **Pick segments** – use the segment dropdowns or accept dependencies pushed
   from other views. Segments aren’t rendered until you click **Render**.
3. **Choose resolution** – either pick a fixed bin size or select **Auto**.
   Auto mode looks at the viewport size and segment lengths, ensuring fewer than
   10 tiles are fetched per render.
4. **Render** – clicking **Render** locks the selection in, fetches any missing
   tiles, and paints the canvas. Changes made after rendering are only applied
   once you click **Render** again, preventing accidental recomputes.

## Features

- **Tile cache** – shared LRU cache keyed by track, segments, resolution, and
  tile coordinates. Keeps re-renders snappy as you zoom/pan.
- **Auto resolution** – calculates the highest resolution that fits the current
  viewport while staying below the tile ceiling. Auto mode refetches tiles only
  after zoom/pan completes.
- **Heatmap export** – downloads a PNG combining the canvas pixels with axis SVG
  to maintain crisp labels.
- **Dependency handling** – listens to cross-view events (e.g., “Generate
  heatmap” from circos) to update its selections without opening the track
  picker manually.

## Tips

- Pair vRISE with Multilift. After mapping loci, send them into a map view to
  immediately visualise interaction patterns in the new reference space.
- Use the **Reset Zoom** button to jump back to the full segment extent instead
  of panning manually.
- When comparing multiple datasets, open multiple map views; the tile cache is
  global, so reusing the same tiles across views stays cheap.
