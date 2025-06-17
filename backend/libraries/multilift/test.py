#!/usr/bin/env python3
"""
multilift_evaluate.py

Sweeps over different synthetic test parameters:
  - number of true intervals
  - insertion length

For each combination, runs N_RUNS repeats of:
  1) Generate random genome + intervals
  2) Insert a block of given length
  3) Run Multilift liftover
  4) Compute TP/FP/FN + precision/recall/F1

Prints results per parameter set.
"""

import os
import random
import tarfile
import tempfile
from io import BytesIO, StringIO

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from main import multilift


# Helper to create FASTA files
def write_fasta(path, seq, id_):
    rec = SeqRecord(Seq(seq), id=id_, description="")
    with open(path, "w") as f:
        SeqIO.write(rec, f, "fasta")


# Generate random genome and intervals
def make_synthetic(seed, length=2000, num_intervals=5):
    random.seed(seed)
    np.random.seed(seed)
    genome = "".join(random.choices("ACGT", k=length))
    starts = sorted(random.sample(range(length), num_intervals))
    intervals = []
    for s in starts:
        e = min(length, s + random.randint(1, 20))
        intervals.append((s, e))
    return genome, intervals


# Produce BED‐formatted bytes
def intervals_to_bytes(intervals):
    txt = "\n".join(f"ref\t{s}\t{e}" for s, e in intervals) + "\n"
    b = txt.encode("utf-8")
    bio = BytesIO(b)
    bio.name = "true_intervals.bed"
    return bio


# Overlap‐based evaluation
def intervals_overlap(a, b):
    return not (a[1] <= b[0] or b[1] <= a[0])


def evaluate(pred, truth):
    TP = FP = FN = 0
    truth_list = [(r.start, r.end) for r in truth.itertuples()]
    pred_list = [(r.start_pos, r.end_pos) for r in pred.itertuples()]
    for p in pred_list:
        if any(intervals_overlap(p, t) for t in truth_list):
            TP += 1
        else:
            FP += 1
    for t in truth_list:
        if not any(intervals_overlap(t, p) for p in pred_list):
            FN += 1
    prec = TP / (TP + FP) if TP + FP > 0 else 0.0
    rec = TP / (TP + FN) if TP + FN > 0 else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0
    return TP, FP, FN, prec, rec, f1


def main():
    N_RUNS = 100
    genome_length = 2000

    interval_settings = [1, 5, 10]
    ins_lengths = [10, 50, 100]

    for num_intervals in interval_settings:
        for INS_LEN in ins_lengths:
            agg = {"TP": 0, "FP": 0, "FN": 0}
            print(f"\n=== Testing num_intervals={num_intervals}, INS_LEN={INS_LEN} ===")
            for run in range(1, N_RUNS + 1):
                seed = 1000 + run

                # 1) generate ref genome + intervals
                work = tempfile.mkdtemp(prefix=f"ml_eval_{run}_")
                gen, intervals = make_synthetic(
                    seed, length=genome_length, num_intervals=num_intervals
                )
                ref_path = os.path.join(work, "ref.fasta")
                mut_path = os.path.join(work, "mut.fasta")
                write_fasta(ref_path, gen, "ref")

                # insert INS_LEN bases at pseudorandom pos
                pos = seed % len(gen)
                ins = "".join(random.choices("ACGT", k=INS_LEN))
                mut_seq = gen[:pos] + ins + gen[pos:]
                write_fasta(mut_path, mut_seq, "mut")

                # 2) build Multilift inputs
                uifmt = ".tar.gz"
                uiupl = False
                genomes = ["ref", "mut"]
                fasta_files = [ref_path, mut_path]
                groups = ["g1"]
                seqs = {}
                for g, fp in zip(genomes, fasta_files):
                    data = open(fp, "r").read()
                    for rec in SeqIO.parse(StringIO(data), "fasta"):
                        seqs[(g, os.path.basename(fp), rec.id, "g1")] = rec

                # 3) wrap true intervals
                bed_io_ref = intervals_to_bytes(intervals)
                bed_io_mut = intervals_to_bytes(intervals)
                uploaded = [bed_io_ref, bed_io_mut]

                # 4) run Multilift
                archive = multilift(
                    uifmt, uiupl, uploaded, groups, seqs, genomes, "mafft"
                )

                # 5) extract prediction
                tar = tarfile.open(fileobj=BytesIO(archive.getvalue()), mode="r:gz")
                member = tar.getmember("liftover/ref/true_intervals.bed")
                f = tar.extractfile(member)
                f.seek(0)
                pred = pd.read_csv(
                    f,
                    sep="\t",
                    header=None,
                    usecols=[0, 1, 2],
                    names=["chrom", "start_pos", "end_pos"],
                )

                # compute shifted ground truth
                shifted = []
                for s, e in intervals:
                    if e <= pos:
                        shifted.append((s, e))
                    elif s >= pos:
                        shifted.append((s + INS_LEN, e + INS_LEN))
                    else:
                        shifted.append((s, e + INS_LEN))
                truth_df = pd.DataFrame(shifted, columns=["start", "end"])

                # evaluate
                TP, FP, FN, prec, rec, f1 = evaluate(pred, truth_df)
                agg["TP"] += TP
                agg["FP"] += FP
                agg["FN"] += FN

            # aggregate
            TP, FP, FN = agg["TP"], agg["FP"], agg["FN"]
            prec = TP / (TP + FP) if TP + FP else 0.0
            rec = TP / (TP + FN) if TP + FN else 0.0
            f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0
            print(f"-> Precision: {prec:.3f}, Recall: {rec:.3f}, F1: {f1:.3f}")


if __name__ == "__main__":
    main()
