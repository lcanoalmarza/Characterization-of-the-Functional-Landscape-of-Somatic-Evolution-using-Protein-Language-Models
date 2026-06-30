#!/usr/bin/env python
# Laura Cano - 14/03/2026

import sys
import os
from Bio import SeqIO

if len(sys.argv) < 2:
    sys.exit("Uso: python script.py <directorio_fastas>")

mane_fasta_dir = sys.argv[1]
ranges = [250, 500, 750, 1000, 10000]
length_ranges = {r: [] for r in ranges}

# Process files
fasta_files = [f for f in os.listdir(mane_fasta_dir) if f.endswith((".fasta", ".faa"))]

for f in fasta_files:
    file_path = os.path.join(mane_fasta_dir, f)
    
    # We only need the length, so we still use next() for efficiency
    record = next(SeqIO.parse(file_path, "fasta"))
    L = len(record.seq)
    
    # Exclusive assignment
    for r in ranges:
        if L <= r:
            # We append the filename 'f' instead of the fasta header 'record.id'
            length_ranges[r].append(f)
            break

# Save results
print(f"{'Range (aa)':<15} | {'File Count':<10}")
print("-" * 28)

for r in ranges:
    files = length_ranges[r]
    print(f"Below {r:<8} | {len(files):<10}")
    
    if files:
        output_path = os.path.join(mane_fasta_dir, f"proteins_range_{r}.txt")
        with open(output_path, "w") as f_out:
            f_out.write("\n".join(files) + "\n")