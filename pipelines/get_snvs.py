#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
get_snvs.py

Simple Python script to extract SNV definitions from a filtered VCF:

 1. Reads a user-specified VCF (skips header lines starting with ‘#’)  
 2. Parses each SNV record (only single-nucleotide substitutions)  
 3. Collects (chrom, pos, REF, ALT) into a list  
 4. Writes that list as JSON to a user-specified output file
"""
import argparse
import json

def main():
    parser = argparse.ArgumentParser(
        description="Extract SNV definitions from a LoFreq-filtered VCF"
    )
    parser.add_argument(
        '-i', '--input', required=True,
        help="Filtered VCF input file"
    )
    parser.add_argument(
        '-o', '--output', required=True,
        help="Output JSON file for SNV list"
    )
    args = parser.parse_args()

    snvs = []
    with open(args.input) as fin:
        for line in fin:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            chrom = fields[0]
            pos   = int(fields[1])
            ref, alt = fields[3], fields[4]
            if len(ref) == 1 and len(alt) == 1:
                snvs.append((chrom, pos, ref, alt))

    with open(args.output, 'w') as fo:
        json.dump(snvs, fo, indent=2)
    print(f"Wrote {len(snvs)} SNVs to {args.output}")


if __name__ == "__main__":
    main()