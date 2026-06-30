# L. Cano - 20/02/2026
# Saturation mutagenesis


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from argparse import ArgumentParser


def get_args():
    """Return the command-line arguments"""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--fasta', help='fasta file with protein sequence(s)')
    parser.add_argument('--output', help='Output file. All corresponding suffixes for single mutation and truncations will be autommatically added')
    return parser.parse_args()

def saturation_mutagenesis(original_record):
    """Yields SeqRecord objects for every possible single-aa mutation"""
    original_seq = str(original_record.seq)
    original_id = original_record.id
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY' # 20 possible aa

    for i, original_aa in enumerate(original_seq):
        for aa in amino_acids:
            if aa == original_aa:
                continue

            # Mutant sequence and id
            mutant_seq = original_seq[:i] + aa + original_seq[i+1:]
            mutant_id = f"{original_id}_{original_aa}{i+1}{aa}"
            yield SeqRecord(Seq(mutant_seq), id=mutant_id, description='')

def generate_truncations(original_record):
    """Yields progressively longer fragments starting from the N-terminus"""
    original_seq = str(original_record.seq)
    original_id = original_record.id

    for i in range(1, len(original_seq) + 1):
        truncated_seq = original_seq[:i]
        mutant_id = f"{original_id}_1_to_{i}"
        yield SeqRecord(Seq(truncated_seq), id=mutant_id, description='')


def main():
    args = get_args()
    
    # Define your paths as strings first
    sm_path = f"{args.output}_SM.fa"
    trunc_path = f"{args.output}_trunc.fa"
    print(f"Input fasta: {args.fasta}")

    # Read input file and compute single-aa mutations and truncations
    with open(sm_path, "w") as f_sm, open(trunc_path, "w") as f_trunc:
        for record in SeqIO.parse(args.fasta, "fasta"):
            # Pass the file OBJECTS (f_sm, f_trunc) to SeqIO.write
            SeqIO.write(saturation_mutagenesis(record), f_sm, "fasta")
            SeqIO.write(generate_truncations(record), f_trunc, "fasta")

    print(f"Files created: {sm_path} and {trunc_path}")


### Main execution
if __name__ == '__main__':
    main()
