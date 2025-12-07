import os
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from hmmlearn import hmm


def run_differential(
    cond1_files, cond2_files, output, *, n_states=2, cov_type="diag", pvalue_cutoff=0.05
):
    """
    cond1_files: list of file paths or file-like objects for condition 1 replicates.
    cond2_files: list of file paths or file-like objects for condition 2 replicates.
    output: file path or file-like to write resulting CSV of significant regions.
    n_states: number of HMM states (default=2).
    cov_type: covariance type for GaussianHMM ('full','diag','tied','spherical'; default='diag').
    pvalue_cutoff: FDR threshold for significance (default=0.05).
    """

    # --- Helper to load a single bedGraph-like replicate, preserving segment info ---
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
            # file-like
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

    # --- Load all replicates for each condition into DataFrames ---
    dfs1 = [load_rep(f, f"rep{i+1}_c1") for i, f in enumerate(cond1_files)]
    dfs2 = [load_rep(f, f"rep{i+1}_c2") for i, f in enumerate(cond2_files)]

    if not dfs1 or not dfs2:
        raise ValueError("At least one replicate per condition is required.")

    # --- Merge replicates within each condition on ['chrom','start','end'] ---
    def merge_all(dfs):
        merged = dfs[0]
        for df in dfs[1:]:
            merged = pd.merge(merged, df, on=["chrom", "start", "end"], how="outer")
        return merged.sort_values(["chrom", "start", "end"]).reset_index(drop=True)

    m1 = merge_all(dfs1)
    m2 = merge_all(dfs2)

    # --- Inner-merge between conditions on ['chrom','start','end'] ---
    combined = pd.merge(m1, m2, on=["chrom", "start", "end"], how="inner")
    combined = combined.sort_values(["chrom", "start", "end"]).reset_index(drop=True)

    # Identify replicate columns
    rep1_cols = [
        c for c in combined.columns if c.startswith("rep") and c.endswith("_c1")
    ]
    rep2_cols = [
        c for c in combined.columns if c.startswith("rep") and c.endswith("_c2")
    ]

    if not rep1_cols:
        raise ValueError("No replicate columns found for condition 1.")
    if not rep2_cols:
        raise ValueError("No replicate columns found for condition 2.")

    # --- Impute missing values (NaNs) with column median ---
    for col in rep1_cols + rep2_cols:
        median_val = combined[col].median(skipna=True)
        combined[col].fillna(median_val, inplace=True)

    # --- Within-replicate normalization: library-size scaling ---
    col_sums = combined[rep1_cols + rep2_cols].sum(axis=0)
    median_sum = np.median(col_sums.values)
    for col in rep1_cols + rep2_cols:
        sf = col_sums[col]
        if sf <= 0:
            continue
        combined[col] = combined[col] * (median_sum / sf)

    # --- Between-replicate normalization: quantile normalization across all replicates ---
    def quantile_normalize(df_reps):
        """
        Perform quantile normalization on a DataFrame where each column is one replicate.
        Returns a DataFrame of the same shape with normalized values.
        """
        sorted_df = np.sort(df_reps.values, axis=0)
        rank_means = np.mean(sorted_df, axis=1)

        df_norm = pd.DataFrame(
            index=df_reps.index, columns=df_reps.columns, dtype=float
        )
        for col in df_reps.columns:
            sorted_idx = np.argsort(df_reps[col].values)
            df_norm[col].iloc[sorted_idx] = rank_means
        return df_norm

    combined_norm_vals = quantile_normalize(combined[rep1_cols + rep2_cols])
    combined[rep1_cols + rep2_cols] = combined_norm_vals

    # --- Compute per-base Δ-reactivity for each row ---
    combined["mean_c1"] = combined[rep1_cols].mean(axis=1)
    combined["mean_c2"] = combined[rep2_cols].mean(axis=1)
    combined["delta"] = np.abs(combined["mean_c2"] - combined["mean_c1"])

    # Collect results across all segments
    results = []

    # --- For each 'chrom' segment, run HMM and region-level tests separately ---
    for chrom, group in combined.groupby("chrom"):
        # 'group' is a DataFrame slice for this segment
        idxs = group.index.to_numpy()  # global indices in 'combined'
        obs_delta = group["delta"].values.reshape(
            -1, 1
        )  # shape = (n_positions_in_segment, 1)

        # Skip if segment is too short
        if obs_delta.shape[0] < 1:
            continue

        # Fit and predict HMM states on Δ-reactivity for this segment
        model = hmm.GaussianHMM(
            n_components=n_states, covariance_type=cov_type, n_iter=100, random_state=0
        )
        model.fit(obs_delta)
        epsilon = 1e-6
        A = model.transmat_
        A += epsilon  # add a tiny pseudocount to every entry
        A /= A.sum(axis=1, keepdims=True)  # renormalize each row to sum to 1
        model.transmat_ = A
        states = model.predict(obs_delta)

        # Determine which HMM state has the higher mean Δ (i.e., "differential" state)
        state_means = model.means_.flatten()
        print(f"Segment: {chrom}, State means: {state_means}")
        peak_state = int(np.argmax(state_means))

        # Extract contiguous runs where state == peak_state
        def segments_from_states(states_array):
            segs = []
            start = None
            for i, s in enumerate(states_array):
                if s == peak_state and start is None:
                    start = i
                elif s != peak_state and start is not None:
                    segs.append((start, i - 1))
                    start = None
            if start is not None:
                segs.append((start, len(states_array) - 1))
            return segs

        segs = segments_from_states(states)

        # For each called region within this segment, compute t-test and effect size
        for loc_start, loc_end in segs:
            # Map local indices to global 'combined' indices
            global_start_idx = idxs[loc_start]
            global_end_idx = idxs[loc_end]

            region_df = group.iloc[loc_start : loc_end + 1].copy()

            # Compute mean per replicate over the region
            region_means_c1 = region_df[rep1_cols].mean(axis=0).values
            region_means_c2 = region_df[rep2_cols].mean(axis=0).values

            region_means_c1 = region_means_c1[~np.isnan(region_means_c1)]
            region_means_c2 = region_means_c2[~np.isnan(region_means_c2)]

            if len(region_means_c1) < 2 or len(region_means_c2) < 2:
                raw_pval = np.nan
            else:
                _, raw_pval = ttest_ind(
                    region_means_c1, region_means_c2, equal_var=False, nan_policy="omit"
                )

            # Compute effect size (mean over all points in region)
            flat_vals_c1 = region_df[rep1_cols].values.flatten()
            flat_vals_c2 = region_df[rep2_cols].values.flatten()
            flat_vals_c1 = flat_vals_c1[~np.isnan(flat_vals_c1)]
            flat_vals_c2 = flat_vals_c2[~np.isnan(flat_vals_c2)]
            if flat_vals_c1.size and flat_vals_c2.size:
                effect_size = np.mean(flat_vals_c2) - np.mean(flat_vals_c1)
            else:
                effect_size = np.nan

            if (
                np.isnan(raw_pval)
                or np.isnan(effect_size)
                or raw_pval >= pvalue_cutoff
                or global_end_idx - global_start_idx <= 1
            ):
                continue

            results.append(
                {
                    "chrom": chrom,
                    "start_pos": combined.loc[global_start_idx, "start"],
                    "end_pos": combined.loc[global_end_idx, "end"],
                    "raw_pval": raw_pval,
                    "effect_size": effect_size,
                }
            )

    # --- Multiple-testing correction (Benjamini–Hochberg) across all segments ---
    df_results = pd.DataFrame(results)
    raw_pvals = df_results["raw_pval"].values
    valid_mask = ~np.isnan(raw_pvals)
    adj_pvals = np.full_like(raw_pvals, np.nan, dtype=float)
    reject_mask = np.zeros_like(raw_pvals, dtype=bool)

    if valid_mask.sum() > 0:
        reject, corrected, _, _ = multipletests(
            raw_pvals[valid_mask], alpha=pvalue_cutoff, method="fdr_bh"
        )
        adj_pvals[valid_mask] = corrected
        reject_mask[valid_mask] = reject

    df_results["adj_pval"] = -np.log10(adj_pvals)
    df_results["is_significant"] = reject_mask

    # --- Write the output CSV of all regions (columns: chrom, start_pos, end_pos, raw_pval, adj_pval, effect_size) ---
    out_columns = [
        "chrom",
        "start_pos",
        "end_pos",
        "raw_pval",
        "adj_pval",
        "effect_size",
    ]
    # out_columns = ['chrom', 'start_pos', 'end_pos', 'raw_pval', 'effect_size']
    if hasattr(output, "write"):
        df_results.to_csv(output, index=False, columns=out_columns)
        output.seek(0)
    else:
        df_results.to_csv(output, index=False, columns=out_columns)

    return df_results


# ---------------------------
# Example CLI stub (unchanged)
# ---------------------------
if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(
        description="Per-segment HMM differential calling on Δ-reactivity"
    )
    p.add_argument(
        "--cond1", nargs="+", required=True, help="Files for condition 1 replicates"
    )
    p.add_argument(
        "--cond2", nargs="+", required=True, help="Files for condition 2 replicates"
    )
    p.add_argument(
        "--states", type=int, default=2, help="Number of HMM states (default: 2)"
    )
    p.add_argument(
        "--cov",
        choices=["full", "diag", "tied", "spherical"],
        default="diag",
        help="Covariance type for GaussianHMM (default: diag)",
    )
    p.add_argument(
        "--pvalue",
        type=float,
        default=0.05,
        help="FDR threshold for significance (default: 0.05)",
    )
    p.add_argument(
        "--output",
        default="diff_peaks.csv",
        help="Output CSV file for significant regions",
    )
    args = p.parse_args()

    run_differential(
        cond1_files=args.cond1,
        cond2_files=args.cond2,
        output=args.output,
        n_states=args.states,
        cov_type=args.cov,
        pvalue_cutoff=args.pvalue,
    )
