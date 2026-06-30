# Simplify multifasta header (retain only id, remove description)
# Laura Cano - 12/02/2026
#!/usr/bin/env python

import sys
import os
from Bio import SeqIO

mane_fasta_dir = sys.argv[1]

fasta_files = [f for f in os.listdir(mane_fasta_dir) if f.endswith((".fasta", ".faa"))]

for fasta_name in fasta_files:
    file_path = os.path.join(mane_fasta_dir, fasta_name)
    
    records = list(SeqIO.parse(file_path, "fasta"))
    
    for record in records:
        record.description = "" 

    with open(file_path, "w") as f:
        SeqIO.write(records, f, "fasta")