#!/usr/bin/env python3
"""Introduce random SNPs into a reference FASTA to create a divergent haplotype.

Used to simulate a heterozygous/diploid test sample for the fungiflow ploidy
detection path (GenomeScope2/Smudgeplot): haplotype A is the original
reference, haplotype B is this script's output at a controlled SNP rate.
"""
import argparse
import random
import sys

from Bio import SeqIO
from Bio.Seq import Seq

BASES = "ACGT"


def mutate_sequence(seq, rate, rng):
    chars = list(str(seq).upper())
    n_mut = 0
    for i, base in enumerate(chars):
        if base not in BASES:
            continue
        if rng.random() < rate:
            choices = [b for b in BASES if b != base]
            chars[i] = rng.choice(choices)
            n_mut += 1
    return Seq("".join(chars)), n_mut


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("reference", help="input reference FASTA (haplotype A)")
    ap.add_argument("output", help="output FASTA (haplotype B)")
    ap.add_argument("--rate", type=float, default=0.015, help="SNP rate per base (default 0.015 = 1.5%%)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--id-suffix", default="_hapB", help="suffix appended to sequence IDs")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    records = list(SeqIO.parse(args.reference, "fasta"))
    total_len = 0
    total_mut = 0
    with open(args.output, "w") as out:
        for rec in records:
            mseq, n_mut = mutate_sequence(rec.seq, args.rate, rng)
            rec.seq = mseq
            rec.id = rec.id + args.id_suffix
            rec.description = ""
            total_len += len(mseq)
            total_mut += n_mut
            SeqIO.write(rec, out, "fasta")

    observed_rate = total_mut / total_len if total_len else 0
    print(
        f"[mutate_haplotype] {len(records)} sequences, {total_len} bp total, "
        f"{total_mut} SNPs introduced, observed rate {observed_rate:.4%} "
        f"(target {args.rate:.4%}), seed={args.seed}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
