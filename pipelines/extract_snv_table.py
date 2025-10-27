#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_snv_table.py

A Python script to convert a filtered VCF into a 7-column TSV variant table.

Steps:
 1. Opens the input VCF and skips header lines.
 2. Parses each SNV record (only single‐nucleotide substitutions).
 3. Extracts INFO fields:
      - DP (total depth)
      - DP4 (ref/alt forward/reverse counts)
      - AF (alternate allele frequency)
 4. Calculates:
      - REF_coverage = DP4[0] + DP4[1]
      - ALT_coverage = DP4[2] + DP4[3]
 5. Writes out the TSV with columns:
      Site, REF, ALT, REF_coverage, ALT_coverage, Depth, ALT_freq
 6. Prints a summary of how many variants were written.
"""

import sys
import argparse

def parse_info_field(info_str):
    d = {}
    for token in info_str.split(';'):
        if '=' in token:
            key, val = token.split('=', 1)
            d[key] = val
    return d

def main():
    p = argparse.ArgumentParser(
        description="Convert a filtered LoFreq VCF into a TSV variant table"
    )
    p.add_argument(
        '-i', '--input',
        required=True,
        help="Filtered VCF file (input)"
    )
    p.add_argument(
        '-o', '--output',
        required=True,
        help="Output TSV file"
    )
    args = p.parse_args()

    try:
        fin = open(args.input, 'r')
    except FileNotFoundError:
        print(f"Error: file {args.input} not found.", file=sys.stderr)
        sys.exit(1)

    fout = open(args.output, 'w')
    fout.write("Site\tREF\tALT\tREF_coverage\tALT_coverage\tDepth\tALT_freq\n")

    written = 0
    for line in fin:
        if line.startswith('#'):
            continue

        fields = line.rstrip('\n').split('\t')
        chrom, pos, vid, ref, alt, qual, filt, info_str = fields[:8]

        # only single‐nucleotide substitutions
        if len(ref) != 1 or len(alt) != 1:
            continue

        info = parse_info_field(info_str)

        dp = int(info.get('DP', '0'))
        dp4_vals = info.get('DP4')
        if not dp4_vals:
            continue
        dp4_list = dp4_vals.split(',')
        if len(dp4_list) < 4:
            continue
        dp4 = [int(x) for x in dp4_list]
        ref_cov = dp4[0] + dp4[1]
        alt_cov = dp4[2] + dp4[3]

        if 'AF' not in info:
            continue
        af = float(info['AF'])

        site_str = f"{chrom}:{pos}"
        fout.write(f"{site_str}\t{ref}\t{alt}\t{ref_cov}\t{alt_cov}\t{dp}\t{af:.4f}\n")
        written += 1

    fin.close()
    fout.close()

    print(f"Wrote {written} variants to file: {args.output}")

if __name__ == "__main__":
    main()