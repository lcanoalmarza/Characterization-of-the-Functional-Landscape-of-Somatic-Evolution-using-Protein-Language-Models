#!/usr/bin/env python
# L. Cano - 20/02/2026
# Exploration of effects of saturation mutagenesis and truncation.

from argparse import ArgumentParser
import numpy as np # type: ignore
import torch
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FuncFormatter, FixedLocator
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
from scipy.ndimage import gaussian_filter1d
from scipy.stats import zscore, kruskal
from pathlib import Path
import requests, sys


def get_args():
    """Return the command-line arguments"""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--wt_fasta', help='wt sequence used to perform Saturation Mutagenesis')
    parser.add_argument('--wt_embedding', help='wt embedding used to perform Saturation Mutagenesis')
    parser.add_argument('--SM_embeddings', help='npz file with protein sequence(s) derived from Saturation Mutagenesis (SM)')
    parser.add_argument('--trunc_embeddings', help='npz file with protein sequence(s) derived from Tuncation Mutagenesis (trunc)')
    parser.add_argument('--model', help='plm used to generate the embeddings')
    parser.add_argument('--output', help='output file with all posible single point mutation')
    parser.add_argument('--sigma', type=int, help='Standard deviation in gaussian filter. Default = 1')
    parser.add_argument('--plot', action='store_true', help='Plot visual exploration')
    parser.add_argument('--umap', action='store_true', help='Perform dimensional reduction of embeddings')

    return parser.parse_args()

def get_embeddings(file_path):
    """
    Read embeddings from .npz (ProtTrans) or a directory of .pt files (ESM).
    Returns a dict of {seq_id: embedding_numpy_array}
    """
    path = Path(file_path)
    embeddings_dict = {}

    # Case 1: directory
    if path.is_dir():
        for pt_file in path.glob("*.pt"):
            data = torch.load(pt_file, map_location="cpu")
            seq_id = data['label']
                
            # Handle both 'mean_representations' and 'representations' keys
            layer_dict = data.get('mean_representations', data.get('representations'))
            if layer_dict:
                layer_idx = next(iter(layer_dict))
                embeddings_dict[seq_id] = layer_dict[layer_idx].float().numpy() # Convert to numpy 


    # Case 2: single .npz file
    elif path.suffix == '.npz':
        with np.load(path) as data:
            ids = data['seq_ids']
            embs = data['embeddings']
            embeddings_dict = {ids[i]: embs[i] for i in range(len(ids))}

    # Case 3: single .pt file 
    elif path.suffix == '.pt':
        data = torch.load(path, map_location="cpu")
        seq_id = data['label']
        layer_dict = data.get('mean_representations', data.get('representations'))
        if layer_dict:
            layer_idx = next(iter(layer_dict))
            embeddings_dict[seq_id] = layer_dict[layer_idx].numpy()

    # Inform (added by Gemini 3.0)
    if not embeddings_dict:
        print(f"Error: No embeddings found at {file_path}")
    else:
        print(f"Successfully loaded {len(embeddings_dict)} embeddings from {path.name}")

    return embeddings_dict

def umap_reduction(embeddings_dict):
    '''Compute UMAP reduction of embeddings'''
    
    # Convert to array
    protein_ids = list(embeddings_dict.keys())
    embeddings_list = [embeddings_dict[pid] for pid in protein_ids]
    X=np.array(embeddings_list)
    
    # Compute UMAP
    reducer=umap.UMAP(metric='cosine', n_components=20)
    red_embeddings=reducer.fit_transform(X)

    # Convet to dict 
    red_embeddings_dict = {pid: red_embeddings[i] for i, pid in enumerate(protein_ids)}
    return red_embeddings_dict

def reduce_all_embeddings(wt_dict, sm_dict, trunc_dict):
    '''Apply UMAP reduction to all embeddings together'''
    
    all_ids = []
    all_vectors = []

    # Process WT 
    for k, v in wt_dict.items():
        all_ids.append(f"WT|{k}")
        all_vectors.append(v)

    # Process SM
    for k, v in sm_dict.items():
        all_ids.append(f"SM|{k}")
        all_vectors.append(v)

    # Process Truncations
    for k, v in trunc_dict.items():
        all_ids.append(f"TR|{k}")
        all_vectors.append(v)

    # Convert to matrix
    X = np.stack(all_vectors)
    
    # Compute UMAP
    print(f"Running UMAP on matrix of shape {X.shape}...")
    reducer = umap.UMAP(metric='cosine', n_components=50, random_state=1984)
    X_reduced = reducer.fit_transform(X)

    # Transform back to dictionaries
    wt_red, sm_red, tr_red = {}, {}, {}
    
    for i, full_id in enumerate(all_ids):
        prefix, original_id = full_id.split("|", 1)
        vector = X_reduced[i]
        
        if prefix == "WT":
            wt_red[original_id] = vector
        elif prefix == "SM":
            sm_red[original_id] = vector
        elif prefix == "TR":
            tr_red[original_id] = vector
            
    return wt_red, sm_red, tr_red

def cosine_distance(a, b):
    """Calculates the cosine distance between two n-dimensional vectors"""
    # Ensure inputs are numpy arrays
    a = np.array(a)
    b = np.array(b)
    
    # Calculate the dot product
    dot_product = np.dot(a, b)
    
    # Calculate the L2 norms (magnitudes)
    norm_a = np.linalg.norm(a)
    norm_b = np.linalg.norm(b)
    
    # Calculate similarity, then distance
    # We add a tiny epsilon (1e-9) to avoid division by zero (Gemini 3)
    similarity = dot_product / (norm_a * norm_b + 1e-9)
    
    return 1 - similarity

def compute_heatmap_SM(wt_fasta, wt_embeddings_dict, SM_embeddings_dict):
    """Compute heatmap of distances between wt and mutated embeddings derived from SM"""
    # Build aa dict
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    aa_dict={ aa:i for aa, i in zip(amino_acids, range(20))}
    
    # Get wt embedding
    wt_embedding=next(iter(wt_embeddings_dict.values()))
    wt_embedding=np.array(wt_embedding).squeeze()

    # Initialize empty matrix
    dist_matrix=np.full((20, len(wt_fasta)), np.nan)

    # Build matrix
    for seq_id, embedding in SM_embeddings_dict.items():
        # Parse headers
        mutation_info = seq_id.split('_')[-1]
        mutant_aa = mutation_info[-1]
        col_idx = int(mutation_info[1:-1]) - 1
        row_idx = aa_dict[mutant_aa]

        # Compute distance wt vs SM
        d = cosine_distance(wt_embedding, embedding)

        # Populate matrix
        dist_matrix[row_idx, col_idx] = d

    return dist_matrix

def normalize_by_aa(dist_matrix):
    """Standardize distances per Amino Acid (row-wise)"""

    # axis=1 calculates the mean/std for each AMINO ACID across all positions
    mu = np.nanmean(dist_matrix, axis=1, keepdims=True)
    sigma = np.nanstd(dist_matrix, axis=1, keepdims=True)
    # Avoid division by zero
    sigma[sigma == 0] = 1.0
    return (dist_matrix - mu) / sigma

#def compute_heatmap_trunc(wt_fasta, wt_embeddings_dict, trunc_embeddings_dict, n_iter=500):
    """Compute heatmap of distances between wt and embeddings of truncated proteins"""

    # Get wt embedding
    wt_embedding=next(iter(wt_embeddings_dict.values()))
    wt_embedding=np.array(wt_embedding).squeeze()

    # Initialize empty vectors
    baseline_vector=np.full((1, len(wt_fasta)), np.nan)     # Averge d of random truncations
    dist_vector=np.full((1, len(wt_fasta)), np.nan)         # d of real truncations
    
    # Build vector
    print(f"Bootstrapping: {n_iter}")
    for seq_id, embedding in trunc_embeddings_dict.items():
        col_idx = int(seq_id.split('_')[-1]) - 1

        # Randomized embeddings distance
        permutation_distances=np.empty((n_iter, len(list(embedding))))
        for i in range(n_iter):
            shuffled_emb = np.random.permutation(embedding)
            d_perm = cosine_distance(wt_embedding, shuffled_emb)
            permutation_distances[i]=d_perm
        baseline_vector[0, col_idx] = np.mean(permutation_distances)

        # Real dist
        d = cosine_distance(wt_embedding, embedding)
        dist_vector[0,col_idx] = (d-baseline_vector[0, col_idx])

    return dist_vector

#def compute_heatmap_trunc(wt_fasta, wt_embeddings_dict, trunc_embeddings_dict):
    """
    Compute heatmap of distances between WT and truncated proteins.
    Uses the distribution of all truncations as the baseline.
    """
    # 1. Get WT embedding and ensure it's a 1D vector (pooled)
    wt_embedding = next(iter(wt_embeddings_dict.values()))
    wt_embedding = np.array(wt_embedding).squeeze()
    if wt_embedding.ndim > 1:
        wt_embedding = wt_embedding.mean(axis=0)

    # Initialize empty vectors
    seq_len = len(wt_fasta)
    dist_vector = np.full((1, seq_len), np.nan)
    
    # 2. First pass: Calculate all real distances
    # We need the mean of ALL truncations to act as the baseline
    temp_dists = {}
    all_d_values = []
    
    for seq_id, embedding in trunc_embeddings_dict.items():
        # Get the index (e.g., from 'seq_123' to 122)
        col_idx = int(seq_id.split('_')[-1]) - 1
        
        # Ensure mutant embedding is pooled
        mut_emb = np.array(embedding).squeeze()
        if mut_emb.ndim > 1:
            mut_emb = mut_emb.mean(axis=0)
            
        d = cosine_distance(wt_embedding, mut_emb)
        
        temp_dists[col_idx] = d
        all_d_values.append(d)

    # 3. Establish the Baseline
    # This represents the 'expected' shift due to length change alone
    baseline_avg = np.mean(all_d_values)
    
    # 4. Second pass: Build the final vector (Observed - Baseline)
    for col_idx, d in temp_dists.items():
        # Positive values = truncation deviates MORE than average
        # Negative values = truncation deviates LESS than average
        dist_vector[0, col_idx] = d - baseline_avg

    return dist_vector

def compute_heatmap_trunc(wt_fasta, wt_embeddings_dict, trunc_embeddings_dict):
    # Get wt embedding - ensure it's a 1D vector (pooled)
    wt_emb = np.array(next(iter(wt_embeddings_dict.values()))).squeeze()
    if wt_emb.ndim > 1: wt_emb = wt_emb.mean(axis=0)

    # Initialize vectors
    seq_len = len(wt_fasta)
    dist_vector = np.full((1, seq_len), np.nan)
    
    # Compute distance of real truncated sequences
    real_dists = {}
    for seq_id, embedding in trunc_embeddings_dict.items():
        col_idx = int(seq_id.split('_')[-1]) - 1
        
        mut_emb = np.array(embedding).squeeze()
        if mut_emb.ndim > 1: mut_emb = mut_emb.mean(axis=0)
        
        real_dists[col_idx] = cosine_distance(wt_emb, mut_emb)

    # Calculate GLOBAL baseline 
    baseline_val = np.mean(list(real_dists.values()))
    
    # Compute normalized embedding deviaiton
    for col_idx, d in real_dists.items():
        dist_vector[0, col_idx] = d - baseline_val

    return dist_vector

def plot_heatmaps(sm_matrix, trunc_vector, amino_acids, p_name, model, output_file):

    # Create figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 10), 
                                   gridspec_kw={'height_ratios': [4, 1]},
                                   sharex=True)
    
    # --- Panel 1: Saturation Mutagenesis ---
    img1 = ax1.imshow(sm_matrix, aspect='auto', cmap='YlOrRd', interpolation='none')
    ax1.set_yticks(range(len(amino_acids)))
    ax1.set_yticklabels(list(amino_acids))
    ax1.set_ylabel("Mutated Amino Acid")
    ax1.set_title(f"Single Amino Acid Saturation Mutagenesis")
    plt.colorbar(img1, ax=ax1, label='Normalised Cosine Distance by Amino Acid')

    # --- Panel 2: Truncation Distance ---
    img2 = ax2.imshow(trunc_vector, aspect='auto', cmap='YlOrRd', interpolation='none')
    ax2.set_yticks([0])
    ax2.set_yticklabels(['Truncation'])
    ax2.set_xlabel("Residue Position")
    ax2.set_title("Protein Truncation")
    plt.colorbar(img2, ax=ax2, label='Cosine Distance')

    # --- X-Axis Ticks Customization ---
    # Major ticks every 50, Minor ticks every 10
    seq_length = sm_matrix.shape[1]
    major_tick_indices = [0] + list(range(9, seq_length, 10))

    def one_based(x, pos):
        return f'{int(x + 1)}'

    for ax in [ax1, ax2]:
        ax.xaxis.set_major_locator(FixedLocator(major_tick_indices))
        ax.xaxis.set_major_formatter(FuncFormatter(one_based))
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        ax.tick_params(axis='x', which='both', labelsize=10, bottom=True, labelbottom=True)

    
    fig.suptitle(f"Saturation Mutagenesis - {p_name.replace('_wt', '')} - {model}", fontsize=22, fontweight='bold', x=0.45, y=0.96)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(f"{output_file}.png", dpi=300)
    plt.close()

def get_uniprot_features(query_val):
    """Retrieve Uniprot Features in a search-like query"""
    # Setup Parameters
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
    "query": query_val,
    "fields": ",".join([
        "accession",
        "protein_name",
        # Functional
        "ft_binding",    # Binding site
        "ft_act_site",   # Active site
        "ft_site",       # Site
                
        # Domains/Regions
        "ft_domain",     # Domain [FT]
        "ft_region",     # Region
        "ft_repeat",     # Repeat
        "ft_motif",      # Motif
        "ft_zn_fing",    # Zinc finger
        "ft_dna_bind", # DNA binding
    
        
        # PTMs / Processing
        "ft_mod_res",    # Modified residue
        "ft_carbohyd",   # Glycosylation 
        "ft_disulfid",   # Disulfide bond
        "ft_signal",     # Signal peptide
        "ft_lipid"       # Lipidation
    ]),
    "sort": "accession desc",
    "size": "1"
    }
    headers = {"accept": "application/json"}

    # Execute Request
    response = requests.get(base_url, headers=headers, params=params)
    if not response.ok:
        response.raise_for_status()
        print("Uniprot features could not be retrieved")
    data = response.json()
    # Parse response
    features_dict = {}
    if data.get('results'):
        protein_data = data['results'][0]
        if 'features' in protein_data:
            for idx, f in enumerate(protein_data['features']):
                f_type = f['type']
                f_start = f['location']['start']['value']
                f_end = f['location']['end']['value']
                f_desc=f['description']
                features_dict[idx] = [f_type, f_start, f_end, f_desc]
    return features_dict

def process_residue_deviation(dist_matrix, output_file, sigma=1):
    # Data processing
    dist_matrix_fixed = np.nan_to_num(dist_matrix, nan=0.0)
    dist_matrix_z = zscore(dist_matrix_fixed, axis=1)
    sensitivity_1d = np.mean(dist_matrix_z, axis=0)
    
    # Log-space smoothing to preserve peaks
    data_to_smooth = sensitivity_1d - np.min(sensitivity_1d) + 1.0
    smoothed_log = gaussian_filter1d(np.log(data_to_smooth), sigma=sigma)
    smoothed_signal = np.exp(smoothed_log) + np.min(sensitivity_1d) - 1.0

    # Save to .npz file
    np.savez(f"{output_file}_per_res_deviation.npz", 
             raw=sensitivity_1d, 
             smoothed=smoothed_signal)
             
    # Save as a simple 2-column TSV 
    np.savetxt(f"{output_file}_per_res_deviation.tsv", 
               np.column_stack((sensitivity_1d, smoothed_signal)), 
               delimiter='\t', header="raw_deviation\tsmoothed_deviation")

    return sensitivity_1d, smoothed_signal

def generate_residue_plot(sensitivity_1d, smoothed_signal, p_name, model, features, output_file, sigma=1):
    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1.5], hspace=0.1)
    ax_main = fig.add_subplot(gs[0])
    ax_track = fig.add_subplot(gs[1], sharex=ax_main)

    color_map = {
        'Domain': '#e67e22', 'Region': "#cbd8e7", 'Motif': "#e73c8c",
        'Zinc finger': "#3cff00", 'DNA binding': "#ff0000",
        'Binding site': "#057434", 'Site': '#95a5a6', 'Modified residue': '#e91e63'
    }

    # Plot deviation signals
    ax_main.plot(sensitivity_1d, alpha=0.15, color='royalblue', label='Raw Deviation')
    ax_main.plot(smoothed_signal, color='darkorange', linewidth=2, label=f'Smoothed (σ={sigma})')
    ax_main.set_ylabel("Avg. Normalised Deviation", fontsize=12)
    ax_main.grid(axis='y', linestyle='--', alpha=0.3)

    # Feature Track Configuration
    feature_types = ['Domain', 'Region', 'Motif', 'Zinc finger', 'DNA binding', 'Binding site']
    occupied_layers = {f_type: [] for f_type in feature_types}
    seen_labels = set()

    if features:
        point_count = 0
        sorted_features = sorted(features.values(), key=lambda x: (x[0], x[1]))

        for f_type, start, end, f_desc in sorted_features:
            color = color_map.get(f_type, "#6fc0fa")
            current_label = f_type if f_type not in seen_labels else None
            if current_label: seen_labels.add(f_type)

            # --- Range Features (Staggered Bars) ---
            if (end - start) > 0 and f_type in occupied_layers:
                layers = occupied_layers[f_type]
                assigned_layer = -1
                for i, layer_spans in enumerate(layers):
                    if not any(max(start-2, s) <= min(end+2, e) for s, e in layer_spans):
                        assigned_layer = i
                        layer_spans.append((start, end))
                        break
                
                if assigned_layer == -1:
                    assigned_layer = len(layers)
                    layers.append([(start, end)])
                
                y_pos = feature_types.index(f_type) * 3 + assigned_layer * 0.8
                ax_track.broken_barh([(start-1, end-start+1)], (y_pos - 0.35, 0.7), 
                                     facecolors=color, alpha=0.8, edgecolor='black', lw=0.5)
                
                if (end - start) > 5:
                    ax_track.text((start + end)/2, y_pos, f_desc[:20], ha='center', 
                                  va='center', fontsize=6, weight='bold', clip_on=True)

            # --- Point Features (Scatter + Annotate) ---
            elif (end - start) == 0:
                pos = start - 1
                if 0 <= pos < len(smoothed_signal):
                    ax_main.scatter(pos, smoothed_signal[pos], color=color, s=40, 
                                    edgecolors='black', zorder=5, label=current_label)
                    offset_y = 10 + (point_count % 3 * 35)
                    ax_main.annotate(f_desc[:35], xy=(pos, smoothed_signal[pos]), 
                                     xytext=(0, offset_y), textcoords='offset points',
                                     fontsize=7, rotation=90, ha='center',
                                     arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
                    point_count += 1
            
            # Legend proxy for range bars
            elif current_label and (end - start) > 0:
                ax_main.scatter([], [], color=color, label=current_label, edgecolor='black', s=40)

    # Styling and Ticks
    ax_track.set_yticks([i * 3 + 0.4 for i in range(len(feature_types))])
    ax_track.set_yticklabels(feature_types, fontsize=9)
    ax_track.set_ylim(-1, len(feature_types) * 3)
    ax_main.set_title(f"{p_name} | Per-Residue Impact | {model}", fontsize=15, pad=20)
    
    seq_len = len(sensitivity_1d)
    ax_track.xaxis.set_major_locator(FixedLocator(np.arange(0, seq_len, 10)))
    ax_track.xaxis.set_major_formatter(FuncFormatter(lambda x, p: f'{int(x+1)}'))
    ax_track.tick_params(axis='x', rotation=90, labelsize=8)
    ax_track.set_xlabel("Residue Position", fontsize=12)
    
    plt.setp(ax_main.get_xticklabels(), visible=False)
    ax_main.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=9)
    
    plt.tight_layout()
    plt.savefig(f"{output_file}_staggered_map.svg", format='svg', bbox_inches='tight')
    plt.show()

def analyze_region_deviation(dist_matrix, features, n_iter, p_name, model, output_file, alpha=0.05):
    '''
    Statistical analysis to determine if embedding deviation is random or
    due to functional disruption in mutant proteins
    '''
    # Data Preparation 
    dist_matrix_fixed = np.nan_to_num(dist_matrix, nan=0.0)
    dist_matrix_z = zscore(dist_matrix_fixed, axis=1)
    sensitivity_1d = np.mean(dist_matrix_z, axis=0)

    # Sort features by start position
    sorted_features = sorted(features.items(), key=lambda x: x[1][1])

    plot_data = []
    labels = []
    long_labels = []
    is_in_region = np.zeros(len(sensitivity_1d), dtype=bool)

    # Feature Processing
    for f_idx, feature in sorted_features:
        f_type, start, end, f_desc = feature
        if start < end:                
            region_values = sensitivity_1d[start-1 : end]
            if region_values.size > 0:
                plot_data.append(region_values)
                long_labels.append(f_desc)
                short_desc = (f_desc[:35] + '..') if len(f_desc) > 35 else f_desc
                labels.append(f"{short_desc}\n({start}-{end})")
                is_in_region[start-1 : end] = True 

    # Add unannotated regions
    unannotated_data = sensitivity_1d[~is_in_region]
    if unannotated_data.size > 0:
        plot_data.append(unannotated_data)
        labels.append('Non-annotated')
        long_labels.append('Non-annotated')

    # --- Kruskal-Wallis test ---
    stat, p_kw = kruskal(*plot_data)
    with open(f"{output_file}_KW_test.tsv", "w") as KW_out:
        KW_out.write(f"H-statistic\t{stat:.4f}\np-value\t{p_kw:.2e}")

    # --- Bootstrapping ---
    num_to_test = len(plot_data) - 1 if labels[-1] == 'Non-annotated' else len(plot_data)
    bootstrap_results = [] # Store p-values for plotting

    with open(f"{output_file}_region_bootsrap.tsv", "w") as bootst_out:
        bootst_out.write("Annotation\tbootsrap-value\tstatus\n")
        
        for i in range(num_to_test):
            region_values = plot_data[i]
            obs_mean = np.mean(region_values)
            
            count_higher = 0
            for _ in range(n_iter):
                random_sample = np.random.choice(sensitivity_1d, size=len(region_values), replace=True)
                if np.mean(random_sample) >= obs_mean:
                    count_higher += 1
            p_boot = count_higher / n_iter
            bootstrap_results.append(p_boot)
            
            # Status Logic
            if p_boot < alpha: status = "*\tUP"
            elif p_boot > (1 - alpha): status = "*\tDOWN"
            else: status = "ns"

            bootst_out.write(f"{long_labels[i]}\t{p_boot:.2e}\t{status}\n")

    # Return everything needed for plotting
    return plot_data, labels, bootstrap_results, sensitivity_1d


#def plot_region_deviation(dist_matrix, p_name, model, features, n_iter, output_file):
    # Data Preparation 
    dist_matrix_fixed = np.nan_to_num(dist_matrix, nan=0.0)
    dist_matrix_z = zscore(dist_matrix_fixed, axis=1)
    sensitivity_1d = np.mean(dist_matrix_z, axis=0)

    # Sort features by start position
    sorted_features = sorted(features.items(), key=lambda x: x[1][1])

    plot_data = []
    labels = []
    long_labels=[]
    is_in_region = np.zeros(len(sensitivity_1d), dtype=bool)

    # Feature Processing
    for f_idx, feature in sorted_features:
        f_type, start, end, f_desc = feature
        if start < end:                
            region_values = sensitivity_1d[start-1 : end]
            if region_values.size > 0:
                plot_data.append(region_values)
                long_labels.append(f_desc)
                short_desc = (f_desc[:35] + '..') if len(f_desc) > 35 else f_desc
                labels.append(f"{short_desc}\n({start}-{end})")
                is_in_region[start-1 : end] = True 

    # Add unannotated regions
    unannotated_data = sensitivity_1d[~is_in_region]
    if unannotated_data.size > 0:
        plot_data.append(unannotated_data)
        labels.append('Non-annotated')

    # Figure Creation 
    fig_height = max(8, len(plot_data) * 0.6)
    fig, ax = plt.subplots(figsize=(12, fig_height))

    
    bplot = ax.boxplot(plot_data, 
                       tick_labels=labels, 
                       patch_artist=True, 
                       vert=False, 
                       showfliers=True,
                       widths=0.6,
                       showmeans=True, 
                       meanline=True) 
    

    # Statistical Analysis (Bootstrap & Kruskal Wallis)
    num_to_test = len(plot_data) - 1 if labels[-1] == 'Non-annotated' else len(plot_data)
    alpha = 0.05 
    
    print(f"\n--- Kruskal-Wallis test for {p_name} ({model}) ---")
    stat, p_kw=kruskal(*plot_data)
    print(f"H-statistic\t{stat:.4f}")
    print(f"p-value\t{p_kw:.2e}")
    with open(f"{output_file}_KW_test.tsv", "w") as KW_out:
        KW_out.write(f"H-statistic\t{stat:.4f}")
        KW_out.write(f"p-value\t{p_kw:.2e}")

    print(f"\n--- Bootstrapping for {p_name} ({model}) ---")
    for i in range(num_to_test):
        region_values = plot_data[i]
        region_len = len(region_values)
        obs_mean = np.mean(region_values)

        # Compute bootstrap replicates
        count_higher = 0
        for _ in range(n_iter):
            random_sample = np.random.choice(sensitivity_1d, size=region_len, replace=True)
            if np.mean(random_sample) >= obs_mean:
                count_higher += 1
        p_boot = count_higher / n_iter
        
        # Save results
        if p_boot < alpha:
            status="*\tUP"
        elif p_boot >  (1-alpha):
            status="*;sDOWN"
        else:
            status=("ns")

        with open(f"{output_file}_region_bootsrap.tsv", "a") as bootst_out:
            bootst_out.write("Annotation\tbootsrap-value\tstatus")
            output_line = f"{long_labels[i].replace('', '')}\t{p_boot:.2e}\t{status}"
            bootst_out.write(output_line + "\n")
        # Mark with asterisk if significant
        if p_boot < alpha or p_boot > 1-alpha:
            x_pos = np.percentile(region_values, 75) + 0.1
            ax.text(x_pos, i + 1, "*", color='red', weight='bold', va='center', fontsize=16)
    

    # Styling 
    colors = plt.cm.RdYlBu_r(np.linspace(0.1, 0.9, len(plot_data)))
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
        patch.set_edgecolor('black')
    ax.set_title(f'{p_name} - Functional Region Impact Analysis\nModel: {model} | Iterations: {n_iter}', pad=20)
    ax.set_xlabel('Average Normalised Deviation (Z-score)')
    ax.axvline(x=0, color='black', linestyle='--', alpha=0.3) 
    ax.grid(axis='x', linestyle=':', alpha=0.5)
    legend_elements = [
        Line2D([0], [0], color='red', lw=2, label='Median'),
        Line2D([0], [0], color='green', lw=2, label='Mean'),
        Line2D([0], [0], marker='o', color='w', label='Outliers', 
               markerfacecolor='gray', markersize=8, alpha=0.6)
    ]

    # Place the legend
    # loc='upper right' or 'lower right' works well; bbox_to_anchor moves it outside
    ax.legend(handles=legend_elements, loc='best', frameon=True, facecolor='white')
    plt.tight_layout()
    plt.savefig(f"{output_file}_region_boxplot.png", dpi=300)
    plt.show()

def generate_region_plot(plot_data, labels, bootstrap_results, p_name, model, n_iter, output_file, alpha=0.05):
    '''
    Boxplot to visually show results from the statistical analysis
    '''
    fig_height = max(8, len(plot_data) * 0.6)
    fig, ax = plt.subplots(figsize=(12, fig_height))

    bplot = ax.boxplot(plot_data, 
                        tick_labels=labels, 
                        patch_artist=True, 
                        vert=False, 
                        showfliers=True,
                        widths=0.6,
                        showmeans=True, 
                        meanline=True) 

    # Mark significance based on pre-calculated bootstrap results
    for i, p_boot in enumerate(bootstrap_results):
        if p_boot < alpha or p_boot > (1 - alpha):
            region_values = plot_data[i]
            x_pos = np.percentile(region_values, 75) + 0.1
            ax.text(x_pos, i + 1, "*", color='red', weight='bold', va='center', fontsize=16)

    # Styling
    colors = plt.cm.RdYlBu_r(np.linspace(0.1, 0.9, len(plot_data)))
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
        patch.set_edgecolor('black')
        
    plt.setp(bplot['medians'], color='red', linewidth=2)
    plt.setp(bplot['means'], color='green', linewidth=2)

    ax.set_title(f'{p_name} - Functional Region Impact\nModel: {model} | Iterations: {n_iter}', pad=20)
    ax.set_xlabel('Average Normalised Deviation (Z-score)')
    ax.axvline(x=0, color='black', linestyle='--', alpha=0.3) 
    ax.grid(axis='x', linestyle=':', alpha=0.5)

    legend_elements = [
        Line2D([0], [0], color='red', lw=2, label='Median'),
        Line2D([0], [0], color='green', lw=2, label='Mean'),
        Line2D([0], [0], marker='o', color='w', label='Outliers', markerfacecolor='gray', markersize=8)
    ]
    ax.legend(handles=legend_elements, loc='best')
    
    plt.tight_layout()
    plt.savefig(f"{output_file}_region_boxplot.png", dpi=300)
    plt.show()

def main():
    # Parse arguments
    args = get_args()
    wt_fasta = args.wt_fasta
    wt_embedding = args.wt_embedding
    SM_embeddings = args.SM_embeddings
    trunc_embeddings = args.trunc_embeddings
    model = args.model
    sigma = args.sigma
    output = args.output
    do_plot = args.plot
    do_umap = args.umap

    # Inform about input files and parameters
    print("-" * 50)
    print(f"{'MUTATIONAL ANALYSIS PIPELINE':^50}")
    print("-" * 50)
    print(f"Wild-Type FASTA:       {wt_fasta}")
    print(f"Wild-Type Embedding:   {wt_embedding}")
    print(f"SM Embeddings:         {SM_embeddings}")
    print(f"Truncation Embeddings: {trunc_embeddings}")
    print("-" * 30)
    print(f"Model used:            {model}")
    print(f"UMAP:                  {do_umap}")
    print(f"Smoothing Sigma:       {sigma}")
    print(f"Output prefix:         {output}")
    print("-" * 50)
    print("\nStarting analysis...\n")

    # Get protein name
    p_name = wt_fasta.split('/')[-1].split('.')[0]

    # Read wt sequence & embeddings
    with open(wt_fasta) as seq:
        wt_sequence = SeqIO.read(seq, 'fasta')
    wt_embeddings_dict = get_embeddings(wt_embedding)
    SM_embeddings_dict = get_embeddings(SM_embeddings)
    trunc_embeddings_dict = get_embeddings(trunc_embeddings)

    # Compute UMAP dimensional reduction
    if do_umap:
        print("Performing UMAP dimensional reduction")
        wt_embeddings_dict, SM_embeddings_dict, trunc_embeddings_dict = reduce_all_embeddings(
            wt_embeddings_dict, SM_embeddings_dict, trunc_embeddings_dict
        )

    # Get uniprot features
    with open(wt_fasta) as f:
        uniprot_search = f.readline().split(".")[0].replace(">", "")
        features_dict = get_uniprot_features(uniprot_search)

    # Compute mutant vs wt distances
    print("Computing embeddings distance respect to wt")
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    SM_dist = compute_heatmap_SM(wt_sequence, wt_embeddings_dict, SM_embeddings_dict)
    norm_SM_dist = normalize_by_aa(SM_dist)
    trunc_dist = compute_heatmap_trunc(wt_sequence, wt_embeddings_dict, trunc_embeddings_dict)
    sens_1d, smooth_sig = process_residue_deviation(norm_SM_dist, output, sigma=sigma)

    # Statistical analysis
    plot_data, labels, bootstrap_results, sensitivity_1d = analyze_region_deviation(
        norm_SM_dist, features_dict, 5000, p_name, model, output, alpha=0.05)
    
    

    # Build plots
    if do_plot:
        plot_heatmaps(norm_SM_dist, trunc_dist, amino_acids, p_name, model, output)
        generate_residue_plot(sens_1d, smooth_sig, p_name, model, features_dict, output, sigma=sigma)
        generate_region_plot(plot_data, labels, bootstrap_results, p_name, model, 5000, output, alpha=0.05)

    print(f"Analysis complete. Results saved with prefix: {output}")
    return



if __name__ == "__main__":
    main()
