#!/usr/bin/env python3
# L. Cano - 26/03/2026
# Compare pLM embedding sensitivity vs AlphaMissense Pathogenicity

import os
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def parse_arguments():
    parser = argparse.ArgumentParser(description="Compare pLMs correlation with AlphaMissense.")
    parser.add_argument("-i", "--input_dir", required=True, help="Directory containing protein folders.")
    parser.add_argument("-f", "--feature_file", required=True, help="TSV file with protein features.")
    parser.add_argument("-a", "--alpha", type=float, default=0.01, help="Significance threshold (default: 0.01).")
    parser.add_argument("-o", "--output", default="alphamissense_comparison.svg", help="Output plot filename.")
    return parser.parse_args()

def main():
    args = parse_arguments()

    # 1. Initialization
    protein_dirs = [f.path for f in os.scandir(args.input_dir) if f.is_dir()]
    proteins = [f.name for f in os.scandir(args.input_dir) if f.is_dir()]
    
    # You can expand this list as you process more models
    model_names = ['ProtT5', 'ProstT5', 'ESM2', 'ESMv1', 'ESMC'] 
    model_colors = ["#EBC033", "#9DB469", "#52A695", "#267FA5", "#08306B"]

    # 2. Read protein types
    protein_type_dict = {}
    with open(args.feature_file) as f_in:
        for line in f_in.readlines()[1:]:
            parts = line.split('\t')
            if len(parts) > 6:
                p_id = parts[0]
                p_type = parts[6].strip() if parts[6].strip() != '' else 'Non-classified'
                protein_type_dict[p_id] = p_type

    # 3. Process data
    comparison_results = []
    for protein, p_dir in zip(proteins, protein_dirs):
        prot_type = protein_type_dict.get(protein, 'Non-classified')
        for model in model_names:
            # Matches the output path format from your previous step
            file_path = os.path.join(p_dir, f"{protein}_{model}_alphamissense_spearman.tsv")
            
            if os.path.exists(file_path):
                try:
                    df_tmp = pd.read_csv(file_path, sep='\t')
                    if not df_tmp.empty:
                        comparison_results.append({
                            'protein': protein,
                            'protein_type': prot_type,
                            'model': model,
                            'rho': float(df_tmp['Rho'].iloc[0]),
                            'p_value': float(df_tmp['p-value'].iloc[0]),
                            'n_residues': int(df_tmp['n_residues'].iloc[0])
                        })
                except Exception as e:
                    print(f"Error reading {file_path}: {e}")

    if not comparison_results:
        print("Error: No data found.")
        return

    df = pd.DataFrame(comparison_results)
    df['is_significant'] = df['p_value'] < args.alpha
    df['model'] = pd.Categorical(df['model'], categories=model_names, ordered=True)

    # 4. Plotting Setup
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial'],
        'pdf.fonttype': 42,         # Standard for academic journals
        'svg.fonttype': 'none',     # CRITICAL: Keeps text as editable text in SVG
        'axes.labelweight': 'bold',
        'axes.linewidth': 1.0,
        'xtick.direction': 'out',
        'ytick.direction': 'out'
    })
    
    sns.set_style("ticks")
    
    protein_types = sorted(df['protein_type'].unique())
    num_types = len(protein_types)
    model_colors = ["#EBC033", "#9DB469", "#52A695", "#267FA5", "#08306B"]

    # Create figure: rows = protein types, cols = 2 (Rho distribution and % Significance)
    fig, axes = plt.subplots(num_types, 2, figsize=(13, 4.5 * num_types), squeeze=False)

    for row_idx, p_type in enumerate(protein_types):
        type_df = df[df['protein_type'] == p_type].copy()
        
        # --- Column A: Distribution of Spearman Rho ---
        sns.violinplot(data=type_df, x='model', y='rho', ax=axes[row_idx, 0], 
                       order=model_names, palette=model_colors, inner=None, alpha=0.6)
        sns.boxplot(data=type_df, x='model', y='rho', ax=axes[row_idx, 0], 
                    order=model_names, width=0.15, color='white', linewidth=1.2, showfliers=False)
        
        axes[row_idx, 0].axhline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
        axes[row_idx, 0].set_ylabel(f"{p_type}\nSpearman ρ", fontsize=10)
        axes[row_idx, 0].set_ylim(-0.3, 1.0) # Standard range for correlation

        # --- Column B: Percentage Significant ---
        perc_df = type_df.groupby('model', observed=True)['is_significant'].mean().reset_index()
        perc_df['percentage'] = perc_df['is_significant'] * 100
        
        sns.barplot(data=perc_df, x='model', y='percentage', ax=axes[row_idx, 1], 
                    order=model_names, palette=model_colors, edgecolor=".2")
        
        axes[row_idx, 1].set_ylabel("Significant Proteins (%)", fontsize=10)
        axes[row_idx, 1].set_ylim(0, 115)

        # Add N counts and % labels on top of bars
        for i, m in enumerate(model_names):
            m_data = type_df[type_df['model'] == m]
            n_total = len(m_data)
            if n_total > 0:
                val = (m_data['is_significant'].mean()) * 100
                n_sig = m_data['is_significant'].sum()
                axes[row_idx, 1].text(i, val + 2, f"{val:.1f}%\n({n_sig}/{n_total})", 
                                     ha='center', va='bottom', fontsize=8, fontweight='bold')

    # Add Panel Labels (a, b) to the top row
    axes[0, 0].set_title('a  Pathogenicity Correlation (Spearman ρ)', loc='left', fontweight='bold', pad=15)
    axes[0, 1].set_title(f'b  Proteins with Stat. Sig. Correlation (α={args.alpha})', loc='left', fontweight='bold', pad=15)

    # Clean up X-axis labels (only show on bottom row)
    for ax in axes.flatten():
        sns.despine(ax=ax)
        ax.set_xlabel('')
    
    plt.tight_layout()

    # Final Save - SVG + PNG for quick preview
    plt.savefig(args.output, format='svg', bbox_inches='tight', transparent=False, facecolor='white')
    
    # Optional: Save a PNG for easy viewing on the cluster
    png_out = args.output.replace('.svg', '.png')
    plt.savefig(png_out, dpi=300, bbox_inches='tight')
    
    plt.close(fig)
    print(f"--- Finished ---")
    print(f"Vector SVG saved: {args.output}")
    print(f"Preview PNG saved: {png_out}")

if __name__ == '__main__':
    main()