#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
consensus_haplotype_freq.py

Enumerates & normalizes all 2^n haplotypes from per‐site allele freqs.

Usage:
  ./consensus_haplotype_freq.py \
    --summary consensus_variant_allele_summary.tsv \
    --output consensus_all_haplotypes.tsv
"""

import argparse
import pandas as pd
import itertools

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--summary', '-s', required=True,
                   help="Per‐site allele summary TSV")
    p.add_argument('--output', '-o', default="consensus_all_haplotypes.tsv",
                   help="Output TSV of full haplotypes")
    args = p.parse_args()

    df = pd.read_csv(args.summary, sep='\t')
    df.set_index('Position', inplace=True)

    allele_probs = {}
    positions = sorted(df.index)
    for pos in positions:
        rec = df.loc[pos]
        allele_probs[pos] = (
            rec['REF_base'], rec['REF_freq'],
            rec['ALT_base'], rec['ALT_freq']
        )

    records = []
    for combo in itertools.product([0,1], repeat=len(positions)):
        hap, prob = [], 1.0
        for i, bit in enumerate(combo):
            pos = positions[i]
            refb, pr, altb, pa = allele_probs[pos]
            if bit == 0:
                hap.append(refb); prob *= pr
            else:
                hap.append(altb); prob *= pa
        if prob>0:
            records.append(("".join(hap), prob))

    hap_df = pd.DataFrame(records, columns=['Haplotype','Frequency'])
    hap_df.sort_values('Frequency', ascending=False, inplace=True)
    hap_df['Frequency'] /= hap_df['Frequency'].sum()
    hap_df.to_csv(args.output, sep='\t', index=False)

    top20 = hap_df.head(20)
    print("Top 20 Haplotypes (Haplotype\tFrequency):")
    for h,f in zip(top20['Haplotype'], top20['Frequency']):
        print(f"{h}\t{f:.6f}")
    print(f"\nTotal haplotypes: {len(hap_df)}")
    print(f"Sum of freqs = {hap_df['Frequency'].sum():.6f}")

if __name__ == "__main__":
    main()