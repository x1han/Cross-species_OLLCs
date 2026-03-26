"""
Cross-Species Cell Type Conservation Analysis using Pseudobulk DESeq2 Results

This module evaluates cell type signature conservation across species using:
1. MetaNeighbor for cell type-level conservation
2. Fisher's Exact Test for marker gene overlap significance
3. Jaccard Index for quantifying conservation at different thresholds
"""

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import pymn
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import fisher_exact
import os

# ============================================================================
# Load Integrated Single-Cell Data
# ============================================================================

adata = sc.read_h5ad('integrated.h5ad')

# ============================================================================
# 1. MetaNeighbor Analysis (Cell Type Level)
# ============================================================================
pymn.variableGenes(adata, study_col='specie')
pymn.MetaNeighborUS(adata, study_col='specie', ct_col='cell_type', fast_version=True)

pymn.plotMetaNeighborUS(adata, figsize=(30, 30), cmap='coolwarm', fontsize=10, show=False)
plt.savefig('MetaNeighbor_AUROC_heatmap.pdf', bbox_inches='tight')
plt.savefig('MetaNeighbor_AUROC_heatmap.png', dpi=600, bbox_inches='tight')
plt.close()

# ============================================================================
# 2. Read Pseudobulk DESeq2 Results
# ============================================================================

# Define path pattern
de_results = {}
for ct in adata.obs['cell_type'].unique():
    de_results[ct] = {}
    for sp in adata.obs['specie'].unique():
        # Read from pseudobulk result files
        file_path = f'lev1_deg/DESeq2/{sp}/{ct}/DESeq2.result.filtered.xls'
        
        if os.path.exists(file_path):
            # Read DESeq2 filtered results (adj.P < 0.05, |log2FC| > 1)
            de_df = pd.read_excel(file_path)
            
            # Extract significantly up-regulated genes as marker genes
            # Note: log2FoldChange > 0 indicates high expression in this cell type vs others
            marker_genes = de_df[de_df['log2FoldChange'] > 0]['gene'].tolist()
            
            de_results[ct][sp] = {
                'markers': set(marker_genes),
                'full_results': de_df
            }
        else:
            de_results[ct][sp] = {'markers': set(), 'full_results': None}

# ============================================================================
# 3. Calculate Jaccard Index and Fisher's Exact Test
# ============================================================================

species_list = list(set([sp for ct in adata.obs['cell_type'].unique() 
                         for sp in de_results[ct].keys() if de_results[ct][sp]['markers']]))
cell_types = adata.obs['cell_type'].unique()

jaccard_results = []
for ct in cell_types:
    for sp1 in species_list:
        for sp2 in species_list:
            if sp1 != sp2:
                genes_sp1 = de_results.get(ct, {}).get(sp1, {}).get('markers', set())
                genes_sp2 = de_results.get(ct, {}).get(sp2, {}).get('markers', set())
                
                if len(genes_sp1) == 0 or len(genes_sp2) == 0:
                    continue
                
                intersection = len(genes_sp1.intersection(genes_sp2))
                union = len(genes_sp1.union(genes_sp2))
                jaccard = intersection / union if union > 0 else 0
                
                # Fisher's exact test
                # 构建列联表
                all_genes_ct = set()
                for sp in de_results[ct].keys():
                    if de_results[ct][sp]['markers']:
                        all_genes_ct.update(de_results[ct][sp]['markers'])
                
                background_size = len(all_genes_ct)
                
                contingency = [[intersection, len(genes_sp1) - intersection],
                               [len(genes_sp2) - intersection, background_size - union]]
                
                try:
                    odds_ratio, p_value = fisher_exact(contingency, alternative='greater')
                except:
                    odds_ratio, p_value = np.nan, np.nan
                
                jaccard_results.append({
                    'cell_type': ct,
                    'species1': sp1,
                    'species2': sp2,
                    'jaccard': jaccard,
                    'odds_ratio': odds_ratio,
                    'p_value': p_value,
                    'shared_genes': intersection,
                    'total_genes_sp1': len(genes_sp1),
                    'total_genes_sp2': len(genes_sp2)
                })

jaccard_df = pd.DataFrame(jaccard_results)
jaccard_df.to_csv('cell_type_conservation.csv', index=False)

# ============================================================================
# 4. Jaccard Index at Different Thresholds (Top N Genes)
# ============================================================================

top_n_list = [200, 500, 1000]
for n in top_n_list:
    results = []
    
    for ct in cell_types:
        for sp1 in species_list:
            for sp2 in species_list:
                if sp1 != sp2:
                    # Get marker genes for each species (sorted by log2FC)
                    markers_sp1 = []
                    markers_sp2 = []
                    
                    if de_results[ct].get(sp1, {}).get('full_results') is not None:
                        df_sp1 = de_results[ct][sp1]['full_results']
                        df_sp1_sorted = df_sp1.sort_values('log2FoldChange', ascending=False)
                        markers_sp1 = set(df_sp1_sorted.head(n)['gene'].tolist())
                    
                    if de_results[ct].get(sp2, {}).get('full_results') is not None:
                        df_sp2 = de_results[ct][sp2]['full_results']
                        df_sp2_sorted = df_sp2.sort_values('log2FoldChange', ascending=False)
                        markers_sp2 = set(df_sp2_sorted.head(n)['gene'].tolist())
                    
                    if len(markers_sp1) == 0 or len(markers_sp2) == 0:
                        continue
                    
                    intersection = len(markers_sp1.intersection(markers_sp2))
                    union = len(markers_sp1.union(markers_sp2))
                    jaccard = intersection / union if union > 0 else 0
                    
                    # Fisher's exact test
                    all_genes_ct = set()
                    for sp in de_results[ct].keys():
                        if de_results[ct][sp]['full_results'] is not None:
                            all_genes_ct.update(de_results[ct][sp]['full_results']['gene'].tolist())
                    
                    background_size = len(all_genes_ct)
                    
                    contingency = [[intersection, len(markers_sp1) - intersection],
                                   [len(markers_sp2) - intersection, background_size - union]]
                    
                    try:
                        odds_ratio, p_value = fisher_exact(contingency, alternative='greater')
                    except:
                        odds_ratio, p_value = np.nan, np.nan
                    
                    results.append({
                        'cell_type': ct,
                        'species1': sp1,
                        'species2': sp2,
                        'jaccard': jaccard,
                        'odds_ratio': odds_ratio,
                        'p_value': p_value,
                        'shared_genes': intersection,
                        'top_n': n
                    })
    
    results_df = pd.DataFrame(results)
    results_df.to_csv(f'cell_type_conservation_top{n}.csv', index=False)
