# Test script to generate mock data for diffStructure
import os
import pandas as pd
import numpy as np
import random
import subprocess
import sys


def load_rep(f, rep_name):
    """
    Load a replicate file in bedGraph-like format with columns:
        chrom<TAB>start<TAB>end<TAB>value
    Returns a DataFrame with columns ['chrom', 'start', 'end', rep_name].
    Filters out rows where value < 0.
    """
    if isinstance(f, str):
        df = pd.read_csv(
            f,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", rep_name],
            usecols=[0, 1, 2, 3],
        )
    else:
        df = pd.read_csv(
            f,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", rep_name],
            usecols=[0, 1, 2, 3],
        )
    # Drop negative values in this replicate
    df = df[df[rep_name] >= 0].copy()
    return df


def generate_mock_data(path: str):
    """
    Input:
        path: str - The path for the original bedgraph data

    Output:
        orig_repls: list - List of original replicates (DataFrames)
        mod_repls: list - List of modified replicates (DataFrames)
        peak_regions: DataFrame - DataFrame containing the pre-selected diff peaks

    Provided an input bedgraph file, the function generates data for two conditions:
        1. Original replicates (orig_repls) with values from the input file.
        2. Modified replicates (mod_repls) with values from the input file, but with
           some values randomly set to simulate more realistic peaks.

    Only one input file is required, and it is used to generate both sets of replicates.
    The function returns two lists of DataFrames, one for original replicates and one for modified replicates.
    """

    # Load the "original" replicate
    original = load_rep(path, "original")
    orig_repls = [original.copy()]
    mod_repls = []

    # --- Create original replicates (4 total) ---
    for i in range(4):
        repl = original.copy().reset_index(drop=True)
        # Rename the "original" column to reflect this replicate
        repl = repl.rename(columns={"original": f"orig_rep{i}"})
        # Randomly perturb some values: with 50% probability, add noise N(0,1)
        repl[f"orig_rep{i}"] = repl[f"orig_rep{i}"].apply(
            # lambda x: x + np.random.normal(0, 1)
            lambda x: x if np.random.rand() < 0.5 else x + np.random.normal(0, 1)
        )
        orig_repls.append(repl)

    # --- Create modified replicates (4 total) ---
    for i in range(4):
        repl = original.copy().reset_index(drop=True)
        # Rename the "original" column to reflect this replicate
        repl = repl.rename(columns={"original": f"mod_rep{i}"})
        # Randomly perturb non-peak positions just like in the original replicates
        repl[f"mod_rep{i}"] = repl[f"mod_rep{i}"].apply(
            lambda x: x if np.random.rand() < 0.5 else x + np.random.normal(0, 1)
            # lambda x: x + np.random.normal(0, 1)
        )
        mod_repls.append(repl)

    # --- Select random "peak" regions ---
    num_regions = len(original)
    peak_regions_list = []
    # Up to 10 peaks, but not more than roughly 10% of total regions
    num_peaks = random.randint(1, min(10, max(1, num_regions // 10)))

    # Choose starting indices for peaks (no overlap enforcement for simplicity)
    chosen_indices = sorted(random.sample(range(num_regions), num_peaks))
    for idx in chosen_indices:
        row = original.iloc[idx]
        chrom = row["chrom"]
        # Pick a random region size between 1 and 10 (or until end of file)
        region_size = random.randint(1, min(10, num_regions - idx))
        end_idx = idx + region_size - 1
        end_idx = min(end_idx, num_regions - 1)
        end_row = original.iloc[end_idx]

        # Draw a peak height: e.g., N(10, 2)
        peak_height = float(np.random.normal(10, 2))

        peak_regions_list.append(
            {
                "chrom": chrom,
                "start": int(row["start"]),
                "end": int(end_row["end"]),
                "peak_height": peak_height,
                "start_idx": idx,
                "end_idx": end_idx,
            }
        )

    peak_regions = pd.DataFrame(peak_regions_list)

    # --- Impose Gaussian‐shaped peaks on the modified replicates ---
    for repl in mod_repls:
        data_col = [c for c in repl.columns if c.startswith("mod_rep")][0]

        for _, region in peak_regions.iterrows():
            chrom = region["chrom"]
            start = region["start"]
            end = region["end"]
            peak_height = region["peak_height"]
            idx0 = region["start_idx"]
            idx1 = region["end_idx"]

            # Mask to find rows falling exactly in this region
            mask = (
                (repl["chrom"] == chrom)
                & (repl["start"] >= start)
                & (repl["end"] <= end)
            )

            # Number of rows in this region
            n = mask.sum()
            if n == 0:
                continue

            # Generate a simple Gaussian shape over these n points
            # Center the peak at the midpoint of [idx0, idx1]
            x = np.arange(n)
            mu = (n - 1) / 2.0
            sigma = n / 4.0  # controls how “wide” the peak is; tweak as desired
            gauss = np.exp(-0.5 * ((x - mu) / sigma) ** 2)

            # Scale so that max(gauss) = peak_height
            gauss_scaled = gauss / gauss.max() * peak_height

            # Add slight per-position noise (e.g. N(0, 0.2))
            noise = np.random.normal(0, 0.2, size=n)
            final_vals = gauss_scaled + noise

            # Now assign these values back into repl[data_col] in true genome order
            # (we assume mask preserves the order)
            repl.loc[mask, data_col] = final_vals

    return orig_repls, mod_repls, peak_regions


def intervals_overlap(a_start, a_end, b_start, b_end):
    return (a_end >= b_start - 80 and a_start <= b_end + 80) or (
        b_end >= a_start - 80 and b_end <= a_end + 80
    )


def evaluate_predictions(pred_df, truth_df):
    # Flatten to lists of intervals per chromosome
    metrics = {"TP": 0, "FP": 0, "FN": 0}
    for chrom in set(truth_df["chrom"]).union(pred_df["chrom"]):
        truth_intervals = truth_df[truth_df["chrom"] == chrom][
            ["start", "end"]
        ].values.tolist()
        pred_intervals = pred_df[pred_df["chrom"] == chrom][
            ["start_pos", "end_pos"]
        ].values.tolist()

        # True positives and false positives
        for p in pred_intervals:
            if any(intervals_overlap(p[0], p[1], t[0], t[1]) for t in truth_intervals):
                metrics["TP"] += 1
            else:
                metrics["FP"] += 1

        # False negatives
        for t in truth_intervals:
            if not any(
                intervals_overlap(t[0], t[1], p[0], p[1]) for p in pred_intervals
            ):
                metrics["FN"] += 1

    # Compute precision, recall, F1
    TP, FP, FN = metrics["TP"], metrics["FP"], metrics["FN"]
    metrics["Precision"] = TP / (TP + FP) if (TP + FP) else 0.0
    metrics["Recall"] = TP / (TP + FN) if (TP + FN) else 0.0
    prec, rec = metrics["Precision"], metrics["Recall"]
    metrics["F1"] = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    return metrics


# if __name__ == "__main__":
#     path = "/Users/mithun/Documents/Imperial/Year 4/FYP/datasets/SHAPE/OC43r11NhighnoUV.xl.bedgraph"
#     orig_repls, mod_repls, peak_regions_df = generate_mock_data(path)

#     output_dir = "/Users/mithun/Documents/Imperial/Year 4/FYP/datasets/SHAPE/mock_data"
#     os.makedirs(output_dir, exist_ok=True)

#     # Save each original replicate as its own bedgraph file
#     for df in orig_repls:
#         rep_name = df.columns[3]
#         out_path = os.path.join(output_dir, f"{rep_name}.bedgraph")
#         df.to_csv(
#             out_path,
#             sep="\t",
#             header=False,
#             index=False,
#             columns=["chrom", "start", "end", rep_name]
#         )

#     # Save each modified replicate as its own bedgraph file
#     for df in mod_repls:
#         rep_name = df.columns[3]
#         out_path = os.path.join(output_dir, f"{rep_name}.bedgraph")
#         df.to_csv(
#             out_path,
#             sep="\t",
#             header=False,
#             index=False,
#             columns=["chrom", "start", "end", rep_name]
#         )

#     # Save the peak regions as a BED file annotation
#     peaks_bed_path = os.path.join(output_dir, "peak_regions.bed")
#     # BED format: chrom, start, end, score (use peak_height as score)
#     peak_regions_df.to_csv(
#         peaks_bed_path,
#         sep="\t",
#         header=False,
#         index=False,
#         columns=["chrom", "start", "end", "peak_height"]
#     )

#     print("BedGraph files for original replicates:")
#     for df in orig_repls:
#         print(f"  - {df.columns[3]}.bedgraph")
#     print("BedGraph files for modified replicates:")
#     for df in mod_repls:
#         print(f"  - {df.columns[3]}.bedgraph")
#     print(f"Peak regions BED written to: {peaks_bed_path}")

#     cond1 = [f"{output_dir}/orig_rep{i}.bedgraph" for i in range(4)]
#     cond2 = [f"{output_dir}/mod_rep{i}.bedgraph"  for i in range(4)]
#     out_csv = f"{output_dir}/diff_peaks.csv"

#     cmd = [
#         sys.executable, "main.py",
#         "--cond1", *cond1,
#         "--cond2", *cond2,
#         "--output", out_csv,
#         "--states", "2",
#         "--cov", "full",
#         "--pvalue", "0.05"
#     ]
#     print("Running:", " ".join(cmd))
#     subprocess.run(cmd, check=True)
#     print("Results written to", out_csv)

#     # 4) load predictions and truth
#     pred = pd.read_csv(out_csv)
#     truth = pd.read_csv(peaks_bed_path, sep="\t", header=None,
#                         names=["chrom","start","end","peak_height"])

#     # 5) evaluate
#     metrics = evaluate_predictions(pred, truth)
#     print("\nEvaluation Metrics:")
#     for k,v in metrics.items():
#         print(f"  {k}: {v:.3f}")

if __name__ == "__main__":
    N_ITER = 10
    states = [2]
    p_vals = [0.01]
    cov = ["diag", "full"]
    total_TP = total_FP = total_FN = 0

    try:
        # print(f"Running with states={state}, p-value={p_val}")
        total_TP = total_FP = total_FN = 0
        for run in range(N_ITER):
            # 1) generate fresh mock data
            path = "/Users/mithun/Documents/Imperial/Year 4/FYP/datasets/SHAPE/OC43r11NhighnoUV.xl.bedgraph"
            orig, mod, truth = generate_mock_data(path)
            output_dir = (
                "/Users/mithun/Documents/Imperial/Year 4/FYP/datasets/SHAPE/mock_data"
            )
            os.makedirs(output_dir, exist_ok=True)

            # write bedGraph files
            for df in orig + mod:
                name = df.columns[3]
                df.to_csv(
                    f"{output_dir}/{name}.bedgraph",
                    sep="\t",
                    header=False,
                    index=False,
                    columns=["chrom", "start", "end", name],
                )
            peaks_bed = f"{output_dir}/peak_regions.bed"
            truth[["chrom", "start", "end", "peak_height"]].to_csv(
                peaks_bed, sep="\t", header=False, index=False
            )

            # 2) call main.py
            cond1 = [f"{output_dir}/orig_rep{i}.bedgraph" for i in range(4)]
            cond2 = [f"{output_dir}/mod_rep{i}.bedgraph" for i in range(4)]
            out_csv = f"{output_dir}/diff_peaks.csv"
            cmd = [
                sys.executable,
                "main.py",
                "--cond1",
                *cond1,
                "--cond2",
                *cond2,
                "--output",
                out_csv,
                "--states",
                "2",
                "--cov",
                "full",
                "--pvalue",
                "0.025",
            ]
            subprocess.run(cmd, check=True)

            # 3) load and evaluate
            pred = pd.read_csv(out_csv)
            truth_df = pd.read_csv(
                peaks_bed,
                sep="\t",
                header=None,
                names=["chrom", "start", "end", "peak_height"],
            )
            metrics = evaluate_predictions(pred, truth_df)

            total_TP += metrics["TP"]
            total_FP += metrics["FP"]
            total_FN += metrics["FN"]
    except Exception as e:
        # 4) compute overall precision, recall, F1
        precision = total_TP / (total_TP + total_FP) if total_TP + total_FP else 0.0
        recall = total_TP / (total_TP + total_FN) if total_TP + total_FN else 0.0
        f1 = (
            2 * precision * recall / (precision + recall) if precision + recall else 0.0
        )

        print(f"Over {N_ITER} runs:")
        print(f"  Total TP: {total_TP}")
        print(f"  Total FP: {total_FP}")
        print(f"  Total FN: {total_FN}")
        print(f"  Overall Precision: {precision:.3f}")
        print(f"  Overall Recall:    {recall:.3f}")
        print(f"  Overall F1 Score:  {f1:.3f}")
    finally:
        # 4) compute overall precision, recall, F1
        precision = total_TP / (total_TP + total_FP) if total_TP + total_FP else 0.0
        recall = total_TP / (total_TP + total_FN) if total_TP + total_FN else 0.0
        f1 = (
            2 * precision * recall / (precision + recall) if precision + recall else 0.0
        )

        print(f"Over {N_ITER} runs:")
        print(f"  Total TP: {total_TP}")
        print(f"  Total FP: {total_FP}")
        print(f"  Total FN: {total_FN}")
        print(f"  Overall Precision: {precision:.3f}")
        print(f"  Overall Recall:    {recall:.3f}")
        print(f"  Overall F1 Score:  {f1:.3f}")
