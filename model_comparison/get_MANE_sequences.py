# Get sequences for the OncoKB list of genes
# Laura Cano - 12/02/2026
#!/usr/bin/env python

import sys
import os
from Bio import SeqIO

oncoKB_path = sys.argv[1]
MANE_path = sys.argv[2]
out_dir = sys.argv[3]

# Create output directory if it doesn't exist
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# Read OncoKB identifiers
with open(oncoKB_path) as f:
    lines = f.readlines()[1:]
    # Clean version numbers from IDs for better matching
    grch38_list = [line.split('\t')[4].split('.')[0] for line in lines]
    hugo_list = [line.split('\t')[0] for line in lines]

# Map MANE records by Gene Symbol and Transcript ID
gene_to_record = {}
fasta_tids = {}

with open(MANE_path) as mane:
    for record in SeqIO.parse(mane, 'fasta'):
        desc = record.description
        
        # Map by Gene Symbol
        if "gene_symbol:" in desc:
            symbol = desc.split("gene_symbol:")[1].split(" ")[0]
            gene_to_record[symbol] = record
            
        # Map by Transcript ID
        if "transcript:" in desc:
            tid = desc.split("transcript:")[1].split(".")[0]
            fasta_tids[tid] = record

# Match and save individual files
saved_count = 0
not_found=list()
for hugo in hugo_list:
    if hugo in gene_to_record:
        record = gene_to_record[hugo]
        # Create filename based on gene symbol
        out_path = os.path.join(out_dir, f"{hugo}.fasta")
        
        with open(out_path, "w") as out_f:
            SeqIO.write(record, out_f, "fasta")
        saved_count += 1
    else:
        not_found.append(record.id)
        
with open(f"{out_dir}/not_found.txt", "w") as f_out:
    for id in not_found:
        f_out.write(f"{id}\n")

print(f"Total genes processed from OncoKB: {len(hugo_list)}")
print(f"Individual FASTA files saved to '{out_dir}': {saved_count}")