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
    def swap_lr_columns(df_subset):
        swap_cols = (
            lambda col: col.replace("L_", "X_").replace("R_", "L_").replace("X_", "R_")
        )
        df_swapped = df_subset.copy()
        df_swapped.columns = [
            swap_cols(c) if c.startswith("L_") or c.startswith("R_") else c
            for c in df_swapped.columns
        ]
        return df_swapped

    correct_start = df[df["L_start"] <= df["R_start"]]
    incorrect_start = df[df["L_start"] > df["R_start"]]
    corrected_df = pd.concat(
        [correct_start, swap_lr_columns(incorrect_start)], ignore_index=True
    )

    correct_seq = corrected_df[corrected_df["L_seqnames"] <= corrected_df["R_seqnames"]]
    incorrect_seq = corrected_df[
        corrected_df["L_seqnames"] > corrected_df["R_seqnames"]
    ]
    reoriented_df = pd.concat(
        [correct_seq, swap_lr_columns(incorrect_seq)], ignore_index=True
    )

    assert (reoriented_df["L_start"] <= reoriented_df["R_start"]).all()
    assert (reoriented_df["L_seqnames"] <= reoriented_df["R_seqnames"]).all()
    assert len(reoriented_df) == len(df)

    return reoriented_df


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
        np.ndarray: Raw or normalised contact matrix, possibly transposed.
    """
    input_flipped = (gene1, gene2) != tuple(sorted([gene1, gene2]))

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
        input_flipped = False
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
        mat = normalise_matrix(mat, len(sub))

    # Flip axes if user passed genes in reverse
    if input_flipped:
        mat = mat.T

    return mat


def generate_contact_map_tile_bedpe(
    bedpe_path,
    gene1,
    gene2,
    fai_path,
    x_start,
    x_end,
    y_start,
    y_end,
    resolution=100,
    normalise=False,
):
    """
    Computes a contact map tile from a BEDPE file for a specific genomic region.

    Parameters:
        bedpe_path (str): Path to BEDPE file.
        gene1 (str): First gene/segment name.
        gene2 (str): Second gene/segment name (same as gene1 if intragenic).
        fai_path (str): Path to the FAI file with gene lengths.
        x_start (int): Start position on gene1 (0-based).
        x_end (int): End position on gene1.
        y_start (int): Start position on gene2 (0-based).
        y_end (int): End position on gene2.
        resolution (int): Bin size in base pairs.
        normalise (bool): Whether to return normalised contact matrix.

    Returns:
        np.ndarray: Contact matrix for the specified tile region.
    """
    # Track if input was flipped
    input_flipped = (gene1, gene2) != tuple(sorted([gene1, gene2]))

    # Always sort genes for consistent processing
    gene1_sorted, gene2_sorted = sorted([gene1, gene2])

    # If genes were flipped, swap the coordinates too
    if input_flipped:
        gene1_sorted, gene2_sorted = gene2_sorted, gene1_sorted
        x_start, y_start = y_start, x_start
        x_end, y_end = y_end, x_end

    # Read FAI
    fai = pd.read_csv(
        fai_path, sep="\t", header=None, names=["gene", "length"], usecols=[0, 1]
    )
    fai_dict = dict(zip(fai["gene"], fai["length"]))

    if gene1_sorted not in fai_dict or gene2_sorted not in fai_dict:
        raise ValueError(
            f"One or both genes not found in FAI: {gene1_sorted}, {gene2_sorted}"
        )

    # Read BEDPE file
    bedpe_cols = ["L_seqnames", "L_start", "L_end", "R_seqnames", "R_start", "R_end"]
    bedpe = pd.read_csv(
        bedpe_path, sep="\t", header=None, usecols=range(6), names=bedpe_cols
    )

    # Filter for intragenic or intergenic
    if gene1_sorted == gene2_sorted:
        sub = bedpe[
            (bedpe["L_seqnames"] == gene1_sorted)
            & (bedpe["R_seqnames"] == gene1_sorted)
        ]
    else:
        # Normalize BEDPE entries so L is always <= R lexicographically
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

        sub = bedpe[
            (bedpe["L_seqnames"] == gene1_sorted)
            & (bedpe["R_seqnames"] == gene2_sorted)
        ]

    if sub.empty:
        print(
            f"No interactions found for tile region {gene1}:{x_start}-{x_end} x {gene2}:{y_start}-{y_end}"
        )
        # Return empty matrix with correct dimensions
        bins_x = int(np.ceil((x_end - x_start) / resolution))
        bins_y = int(np.ceil((y_end - y_start) / resolution))
        return np.zeros((bins_x, bins_y), dtype=int)

    # Filter reads that overlap with the tile region
    # For intragenic case
    if gene1_sorted == gene2_sorted:
        sub = sub[
            (sub["L_start"] < x_end)
            & (sub["L_end"] > x_start)
            & (sub["R_start"] < y_end)
            & (sub["R_end"] > y_start)
        ].copy()
    else:
        # For intergenic, L corresponds to gene1, R to gene2
        sub = sub[
            (sub["L_start"] < x_end)
            & (sub["L_end"] > x_start)
            & (sub["R_start"] < y_end)
            & (sub["R_end"] > y_start)
        ].copy()

    if sub.empty:
        bins_x = int(np.ceil((x_end - x_start) / resolution))
        bins_y = int(np.ceil((y_end - y_start) / resolution))
        return np.zeros((bins_x, bins_y), dtype=int)

    # Generate tile-specific contact map
    mat = get_contact_map_tile_py(
        sub["L_start"].values,
        sub["L_end"].values,
        sub["R_start"].values,
        sub["R_end"].values,
        x_start,
        x_end,
        y_start,
        y_end,
        resolution,
    )

    if normalise:
        mat = normalise_matrix(mat, len(sub))

    # Flip axes if user passed genes in reverse
    if input_flipped:
        mat = mat.T

    return mat


def get_contact_map_tile_py(
    L_start,
    L_end,
    R_start,
    R_end,
    tile_x_start,
    tile_x_end,
    tile_y_start,
    tile_y_end,
    resolution,
):
    """
    Generate a contact map for a specific tile region.

    Parameters:
        L_start, L_end: Arrays of left read start/end positions (1-based from BEDPE)
        R_start, R_end: Arrays of right read start/end positions (1-based from BEDPE)
        tile_x_start, tile_x_end: Genomic coordinates of the tile on x-axis (0-based)
        tile_y_start, tile_y_end: Genomic coordinates of the tile on y-axis (0-based)
        resolution: Bin size in base pairs

    Returns:
        np.ndarray: Contact matrix for the tile region
    """
    # Calculate number of bins in this tile
    bins_x = int(np.ceil((tile_x_end - tile_x_start) / resolution))
    bins_y = int(np.ceil((tile_y_end - tile_y_start) / resolution))

    contact_map = np.zeros((bins_x, bins_y), dtype=int)

    for l1, l2, r1, r2 in zip(L_start, L_end, R_start, R_end):
        # Convert to 0-based
        l1, l2, r1, r2 = l1 - 1, l2 - 1, r1 - 1, r2 - 1

        # Convert genomic coordinates to tile-relative bin coordinates
        # Clamp to tile boundaries
        l1_clamped = max(l1, tile_x_start)
        l2_clamped = min(l2, tile_x_end - 1)
        r1_clamped = max(r1, tile_y_start)
        r2_clamped = min(r2, tile_y_end - 1)

        # Skip if read doesn't overlap tile
        if l1_clamped > l2_clamped or r1_clamped > r2_clamped:
            continue

        # Convert to bin indices relative to tile origin
        bin_x_start = (l1_clamped - tile_x_start) // resolution
        bin_x_end = (l2_clamped - tile_x_start) // resolution
        bin_y_start = (r1_clamped - tile_y_start) // resolution
        bin_y_end = (r2_clamped - tile_y_start) // resolution

        # Fill in the contact map
        for x in range(bin_x_start, min(bin_x_end + 1, bins_x)):
            for y in range(bin_y_start, min(bin_y_end + 1, bins_y)):
                contact_map[x, y] += 1

    return contact_map
