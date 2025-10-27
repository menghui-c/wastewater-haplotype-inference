#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
consensus_summarize_variant_alleles.py

Generates per‐site allele counts and frequencies from a MAFFT‐aligned multi‐FASTA.

Usage:
  ./consensus_summarize_variant_alleles.py \
    --variants-table variants_final_table.tsv \
    --aligned-fasta aligned.fasta \
    --reference-fasta reference.fasta \
    --output consensus_variant_allele_summary.tsv
"""
import argparse
import pandas as pd
from collections import Counter

def read_fasta(filepath):
    seqs = {}
    with open(filepath, 'r') as f:
        cur_id, cur_seq = None, []
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith('>'):
                if cur_id:
                    seqs[cur_id] = ''.join(cur_seq)
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line.upper())
        if cur_id:
            seqs[cur_id] = ''.join(cur_seq)
    return seqs

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--variants-table', '-v', required=True,
                   help="TSV with columns Site,REF,ALT,...")
    p.add_argument('--aligned-fasta', '-a', required=True,
                   help="MAFFT‐aligned FASTA (ref + samples)")
    p.add_argument('--reference-fasta', '-r', required=True,
                   help="Original reference FASTA")
    p.add_argument('--output', '-o', default="consensus_variant_allele_summary.tsv",
                   help="Output TSV path")
    args = p.parse_args()

    df = pd.read_csv(args.variants_table, sep='\t')
    df['Position'] = df['Site'].str.split(':').str[-1].astype(int)
    df['REF_base'] = df['REF']
    df['ALT_base'] = df['ALT']

    seqs = read_fasta(args.aligned_fasta)
    total_seq = len(seqs)
    if total_seq < 2:
        raise RuntimeError("Need ≥1 reference + ≥1 sample in aligned FASTA")
    total_samples = total_seq - 1

    # reference = first sequence in aligned FASTA
    ref_id = list(seqs.keys())[0]
    ref_seq_aligned = seqs[ref_id]

    # load original ref (to compare lengths)
    with open(args.reference_fasta) as f:
        lines = [l.strip() for l in f if not l.startswith('>')]
    ref_seq = "".join(lines).upper()
    if len(ref_seq_aligned) != len(ref_seq):
        raise RuntimeError(
            f"Aligned ref length {len(ref_seq_aligned)} != original {len(ref_seq)}"
        )

    # count alleles per site
    allele_counts = {pos: Counter() for pos in df['Position']}
    for sid, seq in seqs.items():
        for pos in df['Position']:
            allele_counts[pos][ seq[pos-1] ] += 1

    # build summary
    records = []
    for _, row in df.iterrows():
        pos, refb, altb = row['Position'], row['REF_base'], row['ALT_base']
        counts = allele_counts[pos]
        # subtract ref itself
        ref_all = counts.get(refb, 0)
        alt_all = counts.get(altb, 0)
        # sample‐only
        ref_count = max(ref_all - (1 if ref_seq_aligned[pos-1]==refb else 0), 0)
        alt_count = alt_all
        other_count = max(total_samples - (ref_count + alt_count), 0)
        # freqs
        ref_freq, alt_freq = ref_count/total_samples, alt_count/total_samples
        other_freq = other_count/total_samples

        records.append({
            'Position': pos,
            'REF_base': refb,
            'ALT_base': altb,
            'REF_count': ref_count,
            'ALT_count': alt_count,
            'Other_count': other_count,
            'Total_samples': total_samples,
            'REF_freq': round(ref_freq,4),
            'ALT_freq': round(alt_freq,4),
            'Other_freq': round(other_freq,4)
        })

    out_df = pd.DataFrame(records)
    out_df.to_csv(args.output, sep='\t', index=False)
    print(f"[INFO] Written {args.output}; samples = {total_samples}")

    # optionally show via ace_tools
    try:
        import ace_tools as tools
        tools.display_dataframe_to_user(
            name="Allele Counts per SNV Site",
            dataframe=out_df
        )
    except ImportError:
        pass

if __name__ == "__main__":
    main()