#!/usr/bin/env python
# L. Cano - 26/03/2026
# Compare model pathogenicity detection based on comparison with AlphaMissense

import os
import sys
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Compare pLMs embedding sensitivity to single-aa mutations.")
    parser.add_argument("-i", "--input_dir", required=True, help="Directory containing protein folders.")
    parser.add_argument("-f", "--feature_file", required=True, help="TSV file with protein features/classifications.")
    parser.add_argument("-a", "--alpha", type=float, default=0.01, help="Significance threshold (default: 0.01).")
    parser.add_argument("-o", "--output", default="model_comparison.svg", help="Output plot filename.")
    return parser.parse_args()




def main():

    # Parse arguments
    args = parse_arguments()
    input_dir=args.input_dir
    features=args.features_file
    alpha=args.alpha
    output=args.output

    # Get protein list
    protein_dirs = [f.path for f in os.scandir(input_dir) if f.is_dir()]
    proteins = [f.name for f in os.scandir(input_dir) if f.is_dir()]

    print(f"--- Initialization ---")
    model_names = ['ProtT5', 'ProstT5', 'ESM2', 'ESMv1', 'ESMC']

    # Classify proteins per type
    print("Reading feature files...")
    protein_type_dict = {}
    with open(args.feature_file) as f_in:
        for line in f_in.readlines()[1:]:
            parts = line.split('\t')
            if len(parts) > 0:
                p_id = parts[0]
                p_type = parts[6].strip() if len(parts) > 6 and parts[6].strip() != '' else 'Non-classified'
                protein_type_dict[p_id] = p_type

    # Compare models 
    comparison_results = []
    non_comparable_count = 0
    comparable_ids=[]

    print(f"Processing {len(proteins)} proteins...")
    sys.exit()
    alphamissense_spearman=list()

    for protein, p_dir in zip(proteins, protein_dirs):
        prot_type = protein_type_dict.get(protein, 'Non-classified')

        for model in model_names:
            alphamissense_path = os.path.join(p_dir, f"{protein}_{model}_alphamissense_spearman.tsv")
            with open(alphamissense_path) as f_in: 
                for line in f.readlines()[1:]:
                    speearman_result={'protein':protein, 'rho':line.split('\t')[0], 'p_value':line.split('\t')[1], 'n_residues':line.split('\t')[2]}


if __name__== '__main__':
    main()