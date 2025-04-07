import pandas as pd
import numpy as np


def get_contact_map_py(
    L_start, L_end, R_start, R_end, gene_A_size, gene_B_size, resolution
):
    bins_A = gene_A_size // resolution + 1
    bins_B = gene_B_size // resolution + 1
    contact_map = np.zeros((bins_A, bins_B), dtype=int)

    for l1, l2, r1, r2 in zip(L_start, L_end, R_start, R_end):
        l1, l2, r1, r2 = l1 - 1, l2 - 1, r1 - 1, r2 - 1  # Convert to 0-based
        for x in range(l1 // resolution, l2 // resolution + 1):
            for y in range(r1 // resolution, r2 // resolution + 1):
                contact_map[x, y] += 1

    return contact_map


def reorient_hybrids(df):
    incorrect = df[df["L_start"] > df["R_start"]].copy()
    correct = df[df["L_start"] <= df["R_start"]]

    swap_cols = (
        lambda col: col.replace("L_", "X_").replace("R_", "L_").replace("X_", "R_")
    )
    incorrect.columns = [
        swap_cols(c) if c.startswith("L_") or c.startswith("R_") else c
        for c in incorrect.columns
    ]

    df = pd.concat([correct, incorrect], ignore_index=True)
    return df


def normalise_matrix(matrix, hybrids_len):
    """
    Normalise a contact matrix by the number of hybrids.

    Parameters:
        matrix (np.ndarray): Raw contact matrix.
        hybrids_len (int): Number of hybrid reads (used for normalisation).

    Returns:
        np.ndarray: Normalised matrix (same shape).
    """
    norm_matrix = matrix * (1e6 / hybrids_len)
    return norm_matrix


def generate_contact_map(
    hybrids_path, gene1, gene2, fai_path, normalise=False, resolution=100
):
    """
    Computes a contact map for the given gene(s).
    IMPORTANT: Assumes input is a hybrid TSV file with columns:
    - L_seqnames
    - L_start
    - L_end
    - R_seqnames
    - R_start
    - R_end

    Parameters:
        hybrids_path (str): Path to the hybrid TSV file.
        gene1 (str): Name of the first gene.
        gene2 (str): Name of the second gene (can be same as gene1).
        fai_path (str): Path to the FAI file.
        normalise (bool): Whether to return normalised data.

    Returns:
        np.ndarray: Contact map matrix.
    """
    gene1, gene2 = sorted([gene1, gene2], reverse=True)
    fai = pd.read_csv(
        fai_path, sep="\t", header=None, names=["gene", "length"], usecols=[0, 1]
    )
    fai_dict = dict(zip(fai["gene"], fai["length"]))
    hybrids = pd.read_csv(hybrids_path, sep="\t")

    if gene1 not in fai_dict or gene2 not in fai_dict:
        raise ValueError(f"One or both genes not found in FAI: {gene1}, {gene2}")

    if gene1 == gene2:
        genome_size = fai_dict[gene1]
        sub = hybrids[
            (hybrids["type"] == "intragenic") & (hybrids["L_seqnames"] == gene1)
        ]
    else:
        hybrids = reorient_hybrids(hybrids)
        sub = hybrids[
            (hybrids["L_seqnames"] == gene1) & (hybrids["R_seqnames"] == gene2)
        ]
        genome_size = (fai_dict[gene1], fai_dict[gene2])

    if sub.empty:
        print("No hybrids found.")
        return None

    mat = get_contact_map_py(
        sub["L_start"].values,
        sub["L_end"].values,
        sub["R_start"].values,
        sub["R_end"].values,
        genome_size if isinstance(genome_size, int) else genome_size[0],
        genome_size if isinstance(genome_size, int) else genome_size[1],
        resolution,
    )

    if normalise:
        return normalise_matrix(mat, len(sub))
    else:
        return mat


def generate_contact_map_bedpe(
    bedpe_path, gene1, gene2, fai_path, normalise=False, resolution=100
):
    """
    Computes a contact map from a BEDPE file.

    Parameters:
        bedpe_path (str): Path to BEDPE file.
        gene1 (str): First gene name (or segment ID).
        gene2 (str): Second gene name (same as gene1 if intragenic).
        fai_path (str): Path to the FAI file with gene lengths.
        normalise (bool): Whether to return normalised contact matrix.

    Returns:
        np.ndarray: Raw or normalised contact matrix.
    """
    gene1, gene2 = sorted([gene1, gene2])
    fai = pd.read_csv(
        fai_path, sep="\t", header=None, names=["gene", "length"], usecols=[0, 1]
    )
    fai_dict = dict(zip(fai["gene"], fai["length"]))

    if gene1 not in fai_dict or gene2 not in fai_dict:
        raise ValueError(f"One or both genes not found in FAI: {gene1}, {gene2}")

    bedpe_cols = ["L_seqnames", "L_start", "L_end", "R_seqnames", "R_start", "R_end"]
    bedpe = pd.read_csv(
        bedpe_path, sep="\t", header=None, usecols=range(6), names=bedpe_cols
    )

    if gene1 == gene2:
        sub = bedpe[(bedpe["L_seqnames"] == gene1) & (bedpe["R_seqnames"] == gene1)]
        genome_size = fai_dict[gene1]
    else:
        bedpe = bedpe.copy()
        need_swap = bedpe["L_seqnames"] > bedpe["R_seqnames"]
        bedpe.loc[need_swap, ["L_seqnames", "R_seqnames"]] = bedpe.loc[
            need_swap, ["R_seqnames", "L_seqnames"]
        ].values
        bedpe.loc[need_swap, ["L_start", "R_start"]] = bedpe.loc[
            need_swap, ["R_start", "L_start"]
        ].values
        bedpe.loc[need_swap, ["L_end", "R_end"]] = bedpe.loc[
            need_swap, ["R_end", "L_end"]
        ].values

        sub = bedpe[(bedpe["L_seqnames"] == gene1) & (bedpe["R_seqnames"] == gene2)]
        genome_size = (fai_dict[gene1], fai_dict[gene2])

    if sub.empty:
        print("No interactions found.")
        return None

    mat = get_contact_map_py(
        sub["L_start"].values,
        sub["L_end"].values,
        sub["R_start"].values,
        sub["R_end"].values,
        genome_size if isinstance(genome_size, int) else genome_size[0],
        genome_size if isinstance(genome_size, int) else genome_size[1],
        resolution,
    )

    if normalise:
        return normalise_matrix(mat, len(sub))
    else:
        return mat
