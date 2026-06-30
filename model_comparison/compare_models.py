#!/usr/bin/env python3
# Compare pLMs embedding sensitivity to single-aa mutation 
# Updated: 23/03/2026

import os
import sys
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
import numpy as np

def parse_arguments():
    parser = argparse.ArgumentParser(description="Compare pLMs embedding sensitivity to single-aa mutations.")
    parser.add_argument("-i", "--input_dir", required=True, help="Directory containing protein folders.")
    parser.add_argument("-f", "--feature_file", required=True, help="TSV file with protein features/classifications.")
    parser.add_argument("-a", "--alpha", type=float, default=0.01, help="Significance threshold (default: 0.01).")
    parser.add_argument("-o", "--output", default="model_comparison.svg", help="Output plot filename.")
    return parser.parse_args()

def main():
    args = parse_arguments()

    # Get protein list
    protein_dirs = [f.path for f in os.scandir(args.input_dir) if f.is_dir()]
    proteins = [f.name for f in os.scandir(args.input_dir) if f.is_dir()]

    print(f"--- Initialization ---")
    model_names = ['ProtT5', 'ProstT5', 'ESM2', 'ESMv1', 'ESMC', 'VESM_650M']

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

    for protein, p_dir in zip(proteins, protein_dirs):
        prot_type = protein_type_dict.get(protein, 'Non-classified')
        
        for model in model_names:
            kw_path = os.path.join(p_dir, f"{model}_deviation_analysis", f"{model}_KW_test.tsv")
            bootstrap_path = os.path.join(p_dir, f"{model}_deviation_analysis", f"{model}_region_bootsrap.tsv")
            
            if os.path.exists(kw_path) and os.path.exists(bootstrap_path):
                # Parse KW results
                with open(kw_path) as kw_file:
                    lines = kw_file.readlines()
                    h_val = float(lines[0].split('\t')[1].strip())
                    p_val = float(lines[1].split('\t')[1].strip())
            
                # Parse bootstrap results
                with open(bootstrap_path) as bt_file:
                    regions = [line.split('\t')[2].strip() for line in bt_file.readlines()[1:]]
                    func_score = regions.count('*') / len(regions) if regions else 0

                # Save name of protein
                comparison_results.append({
                    'protein': protein, 
                    'protein_type': prot_type, 
                    'model': model, 
                    'H_stat': h_val, 
                    'p_value': p_val, 
                    'functional_score': func_score
                })
            else:
                non_comparable_count += 1

    print(f"Successfully processed entries: {len(comparison_results)}")
    print(f"Missing/Incomplete entries: {non_comparable_count}")


    if not comparison_results:
        print("Error: No data found.")
        return

    # Convert to df
    df = pd.DataFrame.from_records(comparison_results)
    df['is_significant'] = df['p_value'] < args.alpha

    ###### QUANTITATIVE COMPARISON (Global) ######
    pivot_df = df.pivot(index='protein', columns='model', values='H_stat').dropna()
    
    if not pivot_df.empty:
        pairs = list(combinations(model_names, 2))
        p_values = []
        for m1, m2 in pairs:
            _, p = wilcoxon(pivot_df[m1], pivot_df[m2], alternative='two-sided')
            p_values.append(p)
        rejected, corrected_p, _, _ = multipletests(p_values, method='bonferroni')
        print(f"\n--- Statistical Analysis (n={len(pivot_df)} proteins) ---")
        print("Bonferroni correction applied.")

###### PLOTTING (DENSITY PLOT UPDATE) ######
    data_output = args.output.replace(".svg", ".csv")
    df.to_csv(data_output, index=False)

    plt.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial'], 
        'pdf.fonttype': 42, 'axes.labelweight': 'bold'
    })
    sns.set_style("ticks")
    sns.set_context("paper", font_scale=1.1, rc={"axes.linewidth": 1.0})

    protein_types = sorted(df['protein_type'].unique())
    num_types = len(protein_types)
    model_colors = ["#EBC033", "#9DB469", "#52A695", "#267FA5", "#08306B"]

    fig, axes = plt.subplots(num_types, 3, figsize=(18, 4.5 * num_types), squeeze=False)

    for row_idx, p_type in enumerate(protein_types):
        type_df = df[df['protein_type'] == p_type].copy()
        
        # Calculate percentage significance for Plot B
        perc_df = type_df.groupby('model', observed=True)['is_significant'].mean().reset_index()
        perc_df['percentage'] = perc_df['is_significant'] * 100
        perc_df['model'] = pd.Categorical(perc_df['model'], categories=model_names, ordered=True)

        # a. H-stat (Kruskal-Wallis)
        sns.violinplot(data=type_df, x='model', y='H_stat', ax=axes[row_idx, 0], order=model_names, 
                       palette=model_colors, inner=None, alpha=0.7)
        sns.boxplot(data=type_df, x='model', y='H_stat', ax=axes[row_idx, 0], order=model_names, 
                    width=0.12, color='white', linewidth=1.2, showfliers=False)
        
        # b. Significant Proteins (%)
        sns.barplot(data=perc_df, x='model', y='percentage', ax=axes[row_idx, 1], 
                    palette=model_colors, order=model_names, edgecolor=".2")
        axes[row_idx, 1].set_ylim(0, 115)
        
        # c. % Functional Regions Deviated (Density/KDE Plot)
        type_df['func_perc'] = type_df['functional_score'] * 100
        sns.kdeplot(data=type_df, x='func_perc', hue='model', ax=axes[row_idx, 2], 
                    hue_order=model_names, palette=model_colors, fill=True, 
                    common_norm=False, alpha=0.4, linewidth=1.5)
        axes[row_idx, 2].set_xlim(-5, 105)

        # Labels
        axes[row_idx, 0].set_ylabel(f"{p_type}\nH-stat (KW)", fontsize=9)
        axes[row_idx, 1].set_ylabel("Significant Proteins (%)", fontsize=9)
        axes[row_idx, 2].set_ylabel("Density", fontsize=9)
        axes[row_idx, 2].set_xlabel("Deviated Functional Regions (%)", fontsize=9)
        
        # Annotate Plot B
        for i, m in enumerate(model_names):
            row = perc_df[perc_df['model'] == m]
            if not row.empty:
                val = row['percentage'].values[0]
                n_total = len(type_df[type_df['model'] == m])
                n_sig = int(round((val/100) * n_total))
                axes[row_idx, 1].text(i, val + 2, f"{val:.1f}%\n({n_sig}/{n_total})", 
                                      ha='center', va='bottom', fontsize=8, fontweight='bold')

        # Nature-style Panel Titles
        if row_idx == 0:
            axes[row_idx, 0].set_title('a  H-stat Kruskal-Wallis Test', loc='left', fontweight='bold')
            axes[row_idx, 1].set_title(f'b  Significant Effects (α={args.alpha})', loc='left', fontweight='bold')
            axes[row_idx, 2].set_title('c  Functional Regions Deviation Density', loc='left', fontweight='bold')
            # Add legend to the first density plot only to save space
            sns.move_legend(axes[row_idx, 2], "upper right", frameon=False, fontsize=8)
        else:
            # Remove redundant legends for subsequent rows
            if axes[row_idx, 2].get_legend():
                axes[row_idx, 2].get_legend().remove()

    for ax in axes.flatten():
        sns.despine(ax=ax)
        if ax not in axes[:, 2]: # Keep x-label for density plots
            ax.set_xlabel('')

    plt.tight_layout()
    plt.savefig(args.output, format='svg', bbox_inches='tight', dpi=300)
    print(f"Figure saved with density plots: {args.output}")

if __name__ == "__main__":
    main()