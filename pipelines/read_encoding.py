#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
read_encoding.py

General-purpose Python script to:

 1. Load SNV definitions from a JSON file
 2. Open and iterate through all alignments in a BAM
 3. For each read, build a vector of length n (number of SNVs) with:
      “-” if the read does not cover the SNV,
      “0” if it carries the REF allele,
      “1” if it carries the ALT allele
 4. Count the occurrences of each unique pattern across all reads
 5. Write the results to a TSV with columns:
      pattern   read_count

Usage:
    ./read_encoding.py \
      --snv-json snv_list.json \
      --bam input.indelqual.bam \
      --output encoding_counts.tsv
"""

import argparse
import json
import pysam
from collections import Counter

def main():
    parser = argparse.ArgumentParser(
        description="Encode reads at SNV sites into pattern counts"
    )
    parser.add_argument(
        '-s', '--snv-json', required=True,
        help="Input JSON file with SNV definitions"
    )
    parser.add_argument(
        '-b', '--bam', required=True,
        help="Input BAM file (indel-quality recalibrated)"
    )
    parser.add_argument(
        '-o', '--output', required=True,
        help="Output TSV file of pattern counts"
    )
    args = parser.parse_args()

    with open(args.snv_json) as fo:
        snv_list = json.load(fo)

    snvs = [(chrom, pos-1, ref, alt) for chrom, pos, ref, alt in snv_list]

    bam = pysam.AlignmentFile(args.bam, "rb")
    counter = Counter()

    for read in bam.fetch():
        ref2read = {}
        for qpos, rpos in read.get_aligned_pairs(matches_only=True):
            ref2read[rpos] = read.query_sequence[qpos]

        code = []
        for _, pos0, ref, alt in snvs:
            base = ref2read.get(pos0)
            if base is None:
                code.append('-')
            elif base == ref:
                code.append('0')
            elif base == alt:
                code.append('1')
            else:
                code.append('-')
        pattern = ''.join(code)

        if pattern != '-' * len(snvs):
            counter[pattern] += 1

    bam.close()

    with open(args.output, 'w') as fo:
        fo.write("pattern\tread_count\n")
        for pat, cnt in counter.items():
            fo.write(f"{pat}\t{cnt}\n")
    print(f"Wrote {len(counter)} patterns to {args.output}")

if __name__ == "__main__":
    main()
