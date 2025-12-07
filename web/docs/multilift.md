# Multilift Guide

Multilift is TIMGE’s built-in coordinate translation utility. It lets you map
intervals from one reference genome (or assembly version) to another and
quickly visualise the results in other views.

## Overview

| Aspect        | Details |
| ------------- | ------- |
| **Input**     | Source FASTA/2bit file, target FASTA, bed-like list of intervals |
| **Output**    | Mapped intervals in BED, plus optional summary statistics |
| **Best for**  | Comparing orthologous regions, reusing annotations on a new assembly, pre-filtering interesting loci before generating heatmaps |

## Prerequisites

- Upload the source and target FASTA files (FASTA + FAI).
- Ensure the mapping table/liftover chain is accessible to the backend.
- Have a list of segments (BED, CSV, or manual inputs) ready.

## Running Multilift

1. Open the **Multilift** modal from the header button.
2. Select the source reference, target reference, and mapping table.
3. Provide intervals:
   - Paste raw locus strings (e.g. `chr1:1,000-2,000`).
   - Upload a BED file.
   - Pull loci from an existing view’s selection.
4. Click **Run Multilift**. The backend queues a translation job and streams the
   results back as soon as they are ready.

## Interpreting the Output

- **Summary cards** – show how many intervals mapped 1:1, needed splitting, or
  failed. Use this to gauge data quality.
- **Mapped intervals table** – includes source locus, destination locus, and any
  metadata columns carried over from the input BED.
- **Export** – download the translated BED for use in other tools.

## Feeding Results Back

- **Add to map view** – click “Send to Map View” to set the translated loci as
  dependencies. The map view will update its segments to match.
- **Add to circos** – export the mapped BED as a new track file and upload it to
  visualise cross-species relationships.

## FAQs

**Q: What happens if part of a locus fails to map?**  
A: Multilift splits the interval into smaller chunks and reports any failures so
you can inspect them manually.

**Q: Can I use custom mapping tables?**  
A: Yes. Upload them alongside the track files; the backend exposes them under
the same UUID namespace. Update the Multilift form to reference the new file.

**Q: How large can my BED be?**  
A: For interactive use, keep files under ~50k intervals. Larger batches should
either be pre-processed offline or run via a future batch API.
