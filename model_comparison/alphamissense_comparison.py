#!/usr/bin/env python
# L. Cano - Updated 2026
# Explore relation between plm embeddings and AlphaMissense scores

from argparse import ArgumentParser
import numpy as np
import torch
from scipy.stats import spearmanr
from pathlib import Path
from Bio import SeqIO
from collections import defaultdict
import requests

def get_args():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--wt_fasta', required=True)
    parser.add_argument('--wt_embedding', required=True)
    parser.add_argument('--SM_embeddings', required=True)
    parser.add_argument('--alphamissense', required=True)
    parser.add_argument('--model', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--plot', action='store_true')
    return parser.parse_args()

def get_embeddings(file_path):
    path = Path(file_path)
    embeddings_dict = {}
    if path.is_dir():
        for pt_file in path.glob("*.pt"):
            data = torch.load(pt_file, map_location="cpu")
            seq_id = data['label']
            layer_dict = data.get('mean_representations', data.get('representations'))
            if layer_dict:
                layer_idx = next(iter(layer_dict))
                embeddings_dict[seq_id] = layer_dict[layer_idx].to(torch.float32).numpy()
    elif path.suffix == '.npz':
        with np.load(path) as data:
            ids = data['seq_ids']
            embs = data['embeddings']
            embeddings_dict = {ids[i]: embs[i] for i in range(len(ids))}
    return embeddings_dict

def cosine_distance(a, b):
    a, b = np.array(a), np.array(b)
    similarity = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-9)
    return 1 - similarity

def compute_heatmap_SM(wt_seq_record, wt_embeddings_dict, SM_embeddings_dict):
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    aa_dict = {aa: i for i, aa in enumerate(amino_acids)}
    wt_embedding = np.array(next(iter(wt_embeddings_dict.values()))).squeeze()
    dist_matrix = np.full((20, len(wt_seq_record.seq)), np.nan)

    for seq_id, embedding in SM_embeddings_dict.items():
        try:
            mutation_info = seq_id.split('_')[-1]
            mutant_aa = mutation_info[-1]
            col_idx = int(mutation_info[1:-1]) - 1
            if mutant_aa in aa_dict:
                dist_matrix[aa_dict[mutant_aa], col_idx] = cosine_distance(wt_embedding, embedding)
        except Exception:
            continue
    return dist_matrix

def normalize_by_aa(dist_matrix):
    mu = np.nanmean(dist_matrix, axis=1, keepdims=True)
    sigma = np.nanstd(dist_matrix, axis=1, keepdims=True)
    sigma[sigma == 0] = 1.0
    return (dist_matrix - mu) / sigma

def get_uniprot_id(query_val):
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {"query": query_val, "fields": "accession", "size": "1"}
    try:
        response = requests.get(base_url, params=params, timeout=10)
        if response.ok:
            data = response.json()
            if data.get('results'):
                return data['results'][0]['primaryAccession']
    except Exception:
        pass
    return None

def main():
    args = get_args()
    sigma = 1.0 # Default value for display
    
    print(f"\n--- Processing Protein: {Path(args.wt_fasta).stem} ---")
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    # 1. Load Sequences and Embeddings
    wt_seq = SeqIO.read(args.wt_fasta, 'fasta')
    wt_embs = get_embeddings(args.wt_embedding)
    sm_embs = get_embeddings(args.SM_embeddings)

    if not wt_embs or not sm_embs:
        print(f"CRITICAL ERROR: Could not load embeddings for {args.model}")
        return

    # 2. Compute Distances
    dist_matrix = compute_heatmap_SM(wt_seq, wt_embs, sm_embs)
    norm_dist = normalize_by_aa(dist_matrix)
    emb_dev = np.mean(np.nan_to_num(norm_dist, nan=0.0), axis=0)

    # 3. Handle ID mapping and AlphaMissense Loading
    # Use the fasta header as a base ID
    header_id = wt_seq.id.split(".")[0].replace(">", "")
    uniprot_id = get_uniprot_id(header_id)
    
    # Try multiple filename possibilities
    possible_ids = [uniprot_id, header_id, header_id.split('_')[0]]
    ams_file = None
    for pid in possible_ids:
        if pid:
            test_path = Path(args.alphamissense) / f"{pid}.tsv"
            if test_path.exists():
                ams_file = test_path
                break

    if not ams_file:
        print(f"CRITICAL ERROR: AlphaMissense file not found for {header_id} in {args.alphamissense}")
        return

    # 4. Read AlphaMissense Scores
    ams_dict = defaultdict(list)
    with open(ams_file) as f:
        for line in f.readlines()[1:]:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                pos = parts[1][1:-1]
                ams_dict[pos].append(float(parts[2]))

    ams_mean = []
    for i in range(1, len(emb_dev) + 1):
        pos_key = str(i)
        ams_mean.append(np.mean(ams_dict[pos_key]) if pos_key in ams_dict else np.nan)
    ams_mean = np.array(ams_mean)

    # 5. Correlation and Save
    mask = ~np.isnan(ams_mean) & ~np.isnan(emb_dev)
    if np.sum(mask) < 2:
        print("CRITICAL ERROR: No overlapping data points for correlation.")
        return

    rho, p_val = spearmanr(ams_mean[mask], emb_dev[mask])
    
    with open(args.output, "w") as f_out:
        f_out.write("Rho\tp-value\tn_residues\n")
        f_out.write(f"{rho:.3f}\t{p_val:.2e}\t{np.sum(mask)}\n")

    print(f"SUCCESS: Result saved to {args.output}")
    print(f"Spearman Rho: {rho:.3f} | Residues: {np.sum(mask)}")

if __name__ == "__main__":
    main()