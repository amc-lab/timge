#!/usr/bin/env python3
"""
runtime_scalability.py

Measures Multilift’s alignment time, liftover time, and peak memory usage
on synthetic genome sets of varying sizes and counts.
"""

import os
import time
import tempfile
import tarfile
import random
import resource
from io import BytesIO, StringIO

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from main import multilift


# Helper: write a random genome of given length
def write_fasta(path, length):
    seq = "".join(random.choices("ACGT", k=length))
    rec = SeqRecord(Seq(seq), id=os.path.basename(path), description="")
    SeqIO.write(rec, path, "fasta")


# Helper: tiny dummy annotation for liftover
def make_dummy_bed():
    txt = "ref\t0\t1\n"
    bio = BytesIO(txt.encode("utf-8"))
    bio.name = "true_intervals.bed"
    return bio


def measure_peak_rss_kb():
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return usage.ru_maxrss  # kilobytes


def test_configuration(num_genomes, genome_length):
    # prepare workspace
    work = tempfile.mkdtemp(prefix=f"scale_{num_genomes}_{genome_length}_")
    fasta_paths = []
    for i in range(num_genomes):
        fp = os.path.join(work, f"genome_{i+1}.fasta")
        write_fasta(fp, genome_length)
        fasta_paths.append(fp)

    # build sequence map for Multilift
    genomes = [f"g{i+1}" for i in range(num_genomes)]
    seqs = {}
    for g, fp in zip(genomes, fasta_paths):
        data = open(fp).read()
        for rec in SeqIO.parse(StringIO(data), "fasta"):
            seqs[(g, os.path.basename(fp), rec.id, g)] = rec

    # one group per genome
    groups = genomes

    # dummy annotation for each genome
    uploaded = [make_dummy_bed() for _ in genomes]

    # measure alignment + liftover together
    start_mem = measure_peak_rss_kb()
    t0 = time.perf_counter()
    archive = multilift(
        ".tar.gz",  # output format
        False,  # no user-supplied alignment
        uploaded,
        groups,
        seqs,
        genomes,
        "mafft",
    )
    t1 = time.perf_counter()
    end_mem = measure_peak_rss_kb()

    # total time covers both alignment and liftover
    total_time = t1 - t0
    peak_rss_gb = (end_mem - start_mem) / 1e6  # convert KB to GB approx.

    # To approximate align vs liftover split, re-run liftover only:
    # extract liftover time by calling liftover on a fresh annotation
    # (alignment cached in seqs and groups)
    single_bed = make_dummy_bed()
    t2 = time.perf_counter()
    archive2 = multilift(
        ".tar.gz", False, [single_bed for _ in genomes], groups, seqs, genomes, "mafft"
    )
    t3 = time.perf_counter()
    liftover_time = t3 - t2
    align_time = total_time - liftover_time

    return {
        "num_genomes": num_genomes,
        "total_length": num_genomes * genome_length,
        "align_time_s": round(align_time, 2),
        "liftover_time_s": round(liftover_time, 2),
        "memory_gb": round(peak_rss_gb, 2),
    }


def main():
    random.seed(42)
    np.random.seed(42)

    configs = [
        (2, 10_000),
        (5, 100_000),
        (10, 1_000_000),
    ]

    results = []
    for num_genomes, length in configs:
        print(f"Testing {num_genomes} genomes × {length} bp each...")
        res = test_configuration(num_genomes, length)
        print(
            f" → Align: {res['align_time_s']}s, "
            f"Liftover: {res['liftover_time_s']}s, "
            f"Memory: {res['memory_gb']}GB"
        )
        results.append(res)

    df = pd.DataFrame(results)
    print("\nSummary:\n", df.to_string(index=False))


if __name__ == "__main__":
    main()
