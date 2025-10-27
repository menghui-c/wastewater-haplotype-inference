#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compute_haplotype_freqs_allele.py

General-purpose Python script to:
 1) Read patterns & counts from a TSV (`encoding_counts.tsv`)
 2) Read SNV definitions from a JSON (`snv_list.json`)
 3) Automatically detect linkage blocks and compute joint frequencies per block
 4) Calculate marginal allele frequencies at each SNV site
 5) Enumerate all 2^n full haplotypes (REF/ALT strings), compute their frequencies, and sort descending
 6) Write output TSV and confirm the frequencies sum to 1

Usage:
    ./compute_haplotype_freqs_allele.py \
        --counts encoding_counts.tsv \
        --snv-json snv_list.json \
        --output all_hap_allele_freqs.tsv
"""

import argparse
import itertools
import json
import pandas as pd
from collections import defaultdict, deque

def main():
    p = argparse.ArgumentParser(description="Compute full‐genome haplotype frequencies from allele‐encoded reads")
    p.add_argument('-c', '--counts',   required=True, help="TSV file with columns 'pattern' and 'read_count'")
    p.add_argument('-s', '--snv-json', required=True, help="JSON file with SNV list [[chrom, pos, REF, ALT], ...]")
    p.add_argument('-o', '--output',   required=True, help="Output TSV file for haplotype frequencies")
    args = p.parse_args()

    # 1) Read patterns & counts
    df = pd.read_csv(args.counts, sep="\t")
    if df.empty:
        raise RuntimeError(f"No patterns read from {args.counts}")
    hap_counts = dict(zip(df['pattern'], df['read_count']))
    first_pat = df['pattern'].iat[0]
    n = len(first_pat)

    # 2) Read SNV definitions
    snvs = json.load(open(args.snv_json))
    refs = [r for _, _, r, _ in snvs]
    alts = [a for _, _, _, a in snvs]

    # 3) Marginal counts per site
    ref_cnt = [0]*n
    alt_cnt = [0]*n
    for pat, cnt in hap_counts.items():
        for i, c in enumerate(pat):
            if c == '0': ref_cnt[i] += cnt
            elif c == '1': alt_cnt[i] += cnt

    p_ref = []
    p_alt = []
    for i in range(n):
        total = ref_cnt[i] + alt_cnt[i]
        if total > 0:
            p_ref.append(ref_cnt[i]/total)
            p_alt.append(alt_cnt[i]/total)
        else:
            p_ref.append(0.5)
            p_alt.append(0.5)

    # 4) Detect linkage blocks
    graph = {i:set() for i in range(n)}
    for pat in hap_counts:
        covered = [i for i, c in enumerate(pat) if c != '-']
        for i,j in itertools.combinations(covered, 2):
            graph[i].add(j); graph[j].add(i)
    visited = [False]*n
    blocks = []
    for i in range(n):
        if not visited[i]:
            queue = deque([i]); visited[i]=True; comp=[]
            while queue:
                u = queue.popleft(); comp.append(u)
                for v in graph[u]:
                    if not visited[v]:
                        visited[v]=True; queue.append(v)
            blocks.append(comp)

    # 5) Block‐wise frequencies
    block_freqs = []
    for block in blocks:
        cnt = defaultdict(int)
        for pat, c in hap_counts.items():
            sub = ''.join(pat[i] for i in block)
            if '-' not in sub: cnt[sub] += c
        total = sum(cnt.values())
        if total > 0:
            freq = {h: cnt[h]/total for h in cnt}
        else:
            freq = {}
            for bits in itertools.product('01', repeat=len(block)):
                prod = 1.0
                for idx, b in zip(block, bits):
                    prod *= (p_ref[idx] if b=='0' else p_alt[idx])
                freq[''.join(bits)] = prod
            S = sum(freq.values())
            for h in freq: freq[h] /= S
        block_freqs.append((block, freq))

    # 6) Enumerate full haplotypes
    results = []
    for bits in itertools.product('01', repeat=n):
        f = 1.0
        for i, b in enumerate(bits):
            f *= (p_ref[i] if b=='0' else p_alt[i])
        for block, freq_dict in block_freqs:
            key = ''.join(bits[i] for i in block)
            denom = 1.0
            for i in block:
                denom *= (p_ref[i] if bits[i]=='0' else p_alt[i])
            f *= (freq_dict.get(key, 0.0) / denom)
        alleles = ''.join(refs[i] if b=='0' else alts[i] for i,b in enumerate(bits))
        results.append((alleles, f))

    results.sort(key=lambda x: x[1], reverse=True)

    # 7) Write output
    total_f = 0.0
    with open(args.output, 'w') as fo:
        fo.write('haplotype\tfrequency\n')
        for hap, freq in results:
            fo.write(f"{hap}\t{freq:.8e}\n")
            total_f += freq

    print(f"Sum of frequencies = {total_f:.6f}")

if __name__ == "__main__":
    main()