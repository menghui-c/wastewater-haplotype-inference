#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
filter_lofreq_vcf.py

A script to perform secondary SNV filtering on a LoFreq VCF.

Reads a VCF from the first filtering pass, retains only SNVs passing
depth and allele-frequency thresholds, and writes a second-pass VCF.

Usage:
    ./filter_lofreq_vcf.py \
      --input  first_pass.lofreq.raw.vcf \
      --output second_pass.lofreq.filtered.vcf \
      [--min-cov 30] \
      [--min-alt-freq 0.10] \
      [--min-ref-freq 0.10]
"""

import sys
import argparse

def compute_ref_freq(dp4, dp):
    ref_cov = dp4[0] + dp4[1]
    return (ref_cov / dp) if dp > 0 else 0.0

def parse_info_field(info_str):
    d = {}
    for token in info_str.split(';'):
        if '=' in token:
            key, val = token.split('=', 1)
            d[key] = val
    return d

def main():
    parser = argparse.ArgumentParser(
        description="Perform secondary SNV filtering on a LoFreq VCF"
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help="First-pass LoFreq VCF input file"
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help="Second-pass filtered VCF output file"
    )
    parser.add_argument(
        '--min-cov', type=int, default=30,
        help="Minimum total depth (DP) [default: 30]"
    )
    parser.add_argument(
        '--min-alt-freq', type=float, default=0.10,
        help="Minimum ALT allele frequency (AF) [default: 0.10]"
    )
    parser.add_argument(
        '--min-ref-freq', type=float, default=0.10,
        help="Minimum REF allele frequency (computed from DP4) [default: 0.10]"
    )
    args = parser.parse_args()

    try:
        fin = open(args.input, 'r')
    except FileNotFoundError:
        print(f"Error: file {args.input} not found.", file=sys.stderr)
        sys.exit(1)
    fout = open(args.output, 'w')

    total = 0
    kept = 0
    for line in fin:
        if line.startswith('#'):
            fout.write(line)
            continue

        total += 1
        fields = line.rstrip('\n').split('\t')
        chrom, pos, vid, ref, alt, qual, filt, info_str = fields[:8]

        # only SNVs
        if len(ref) != 1 or len(alt) != 1:
            continue

        info = parse_info_field(info_str)

        dp = int(info.get('DP', '0'))
        if dp < args.min_cov:
            continue

        af = float(info.get('AF', 0.0))
        if af < args.min_alt_freq:
            continue

        dp4_vals = info.get('DP4', '').split(',')
        if len(dp4_vals) < 4:
            continue
        dp4 = [int(x) for x in dp4_vals]
        ref_freq = compute_ref_freq(dp4, dp)
        if ref_freq < args.min_ref_freq:
            continue

        fout.write(line)
        kept += 1

    fin.close()
    fout.close()

    print(f"Processed {total} candidate SNVs; retained {kept} passing variants.")
    print(f"Filtered VCF written to: {args.output}")

if __name__ == "__main__":
    main()