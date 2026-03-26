import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

# ============================================================================
# Aging and Disease Risk Gene Analysis (Integrated with Pseudobulk Results)
# ============================================================================

adata = sc.read_h5ad('integrated.h5ad')

# ============================================================================
# 1. Load Risk Gene Lists
# ============================================================================

# Load risk genes from databases
risk_genes_dict = {}

# Aging genes
if os.path.exists('aging_genes.csv'):
    aging_genes = pd.read_csv('aging_genes.csv')['gene'].tolist()
    risk_genes_dict['aging'] = set(aging_genes)

# Disease genes (from DisGeNET or other databases)
disease_files = {
    'AD': 'AD_genes.csv',
    'PD': 'PD_genes.csv', 
    'MS': 'MS_genes.csv',
    'HD': 'HD_genes.csv',
    'SCZ': 'SCZ_genes.csv',
    'ASD': 'ASD_genes.csv',
    'MDD': 'MDD_genes.csv'
}

for disease, file_path in disease_files.items():
    if os.path.exists(file_path):
        disease_genes = pd.read_csv(file_path)['gene'].tolist()
        risk_genes_dict[disease] = set(disease_genes)

# ============================================================================
# 2. Read Aging Pseudobulk Results
# ============================================================================

aging_de_results = {}
cell_types_aging = ['OPC', 'OL']  # Or finer subtypes
species_list = ['human', 'fasci', 'mulatta', 'mouse']

for sp in species_list:
    aging_de_results[sp] = {}
    for ct in cell_types_aging:
        file_path = f'aging/{sp}/{ct}/DESeq2.result.filtered.xls'
        
        if os.path.exists(file_path):
            de_df = pd.read_excel(file_path)
            
            # Extract significantly differentially expressed genes
            deg_up = de_df[de_df['log2FoldChange'] > 0]['gene'].tolist()
            deg_down = de_df[de_df['log2FoldChange'] < 0]['gene'].tolist()
            
            aging_de_results[sp][ct] = {
                'deg_up': set(deg_up),
                'deg_down': set(deg_down),
                'all_deg': set(de_df['gene'].tolist()),
                'full_results': de_df
            }
        else:
            aging_de_results[sp][ct] = {'deg_up': set(), 'deg_down': set(), 'all_deg': set(), 'full_results': None}

# ============================================================================
# 3. Calculate Disease/Aging Risk Module Scores
# ============================================================================
# 3. Calculate Disease/Aging Risk Module Scores
# ============================================================================

for gene_set_name, gene_list in risk_genes_dict.items():
    # Filter genes present in the dataset
    gene_list_filtered = [g for g in gene_list if g in adata.var_names]
    
    if len(gene_list_filtered) < 10:
        continue
    
    # Calculate module scores using score_genes
    try:
        sc.tl.score_genes(adata, gene_list_filtered, score_name=f'{gene_set_name}_score')
    except Exception as e:
        pass

# Save AnnData object with risk scores
adata.write_h5ad('integrated_with_risk_scores.h5ad')

# ============================================================================
# 4. Cross-species Expression Correlation Analysis
# ============================================================================

correlation_results = []
cell_types = adata.obs['cell_type'].unique()

for gene_set in risk_genes_dict.keys():
    if f'{gene_set}_score' not in adata.obs.columns:
        continue
    
    for sp1 in species_list:
        for sp2 in species_list:
            if sp1 == sp2:
                continue
            
            for ct in cell_types:
                # Get risk scores for this cell type in both species
                mask1 = (adata.obs['specie'] == sp1) & (adata.obs['cell_type'] == ct)
                mask2 = (adata.obs['specie'] == sp2) & (adata.obs['cell_type'] == ct)
                
                scores1 = adata.obs.loc[mask1, f'{gene_set}_score']
                scores2 = adata.obs.loc[mask2, f'{gene_set}_score']
                
                if len(scores1) > 0 and len(scores2) > 0:
                    # Calculate correlation of mean expression levels
                    mean_score1 = scores1.mean()
                    mean_score2 = scores2.mean()
                    
                    # Spearman correlation (simplified to single value comparison)
                    # For gene-level correlation, should calculate expression correlation per gene
                    corr, pval = stats.spearmanr([mean_score1], [mean_score2])
                    
                    correlation_results.append({
                        'gene_set': gene_set,
                        'species1': sp1,
                        'species2': sp2,
                        'cell_type': ct,
                        'correlation': corr,
                        'pvalue': pval,
                        'mean_score_sp1': mean_score1,
                        'mean_score_sp2': mean_score2
                    })

correlation_df = pd.DataFrame(correlation_results)
correlation_df.to_csv('disease_correlation.csv', index=False)

# ============================================================================
# 5. Risk Gene Classification (Based on Pseudobulk Results)
# ============================================================================

risk_gene_classification = []

for gene_set_name, gene_list in risk_genes_dict.items():
    for gene in gene_list:
        # Check expression pattern of this gene across different species
        
        classification = {
            'gene': gene,
            'gene_set': gene_set_name,
            'conserved_across_all': False,
            'primate_specific': False,
            'human_specific': False,
            'mouse_specific': False,
            'present_in_species': []
        }
        
        # Check in which species the gene is expressed (based on pseudobulk results)
        present_in = []
        for sp in species_list:
            # Check if expressed in at least one cell type
            for ct in cell_types_aging:
                if aging_de_results.get(sp, {}).get(ct, {}).get('full_results') is not None:
                    df_sp = aging_de_results[sp][ct]['full_results']
                    if gene in df_sp['gene'].values:
                        present_in.append(sp)
                        break
        
        classification['present_in_species'] = list(set(present_in))
        
        # Classification criteria
        if len(set(present_in)) == len(species_list):
            classification['conserved_across_all'] = True
        elif all(sp in present_in for sp in ['human', 'fasci', 'mulatta']) and 'mouse' not in present_in:
            classification['primate_specific'] = True
        elif present_in == ['human']:
            classification['human_specific'] = True
        elif present_in == ['mouse']:
            classification['mouse_specific'] = True
        
        risk_gene_classification.append(classification)

classification_df = pd.DataFrame(risk_gene_classification)
classification_df.to_csv('risk_gene_classification.csv', index=False)

# Summarize classification results
summary_stats = {
    'conserved_across_all': classification_df['conserved_across_all'].sum(),
    'primate_specific': classification_df['primate_specific'].sum(),
    'human_specific': classification_df['human_specific'].sum(),
    'mouse_specific': classification_df['mouse_specific'].sum()
}

# ============================================================================
# 6. Prepare Visualization Data
# ============================================================================

# Correlation heatmap data
heatmap_data = correlation_df.pivot_table(
    values='correlation', 
    index=['gene_set', 'species1'], 
    columns='species2'
)

# Save data for visualization
heatmap_data.to_csv('disease_correlation_heatmap_input.csv')

# Classification pie chart data
classification_summary = []
for gene_set in risk_genes_dict.keys():
    subset = classification_df[classification_df['gene_set'] == gene_set]
    classification_summary.append({
        'gene_set': gene_set,
        'total_genes': len(subset),
        'conserved': subset['conserved_across_all'].sum(),
        'primate_specific': subset['primate_specific'].sum(),
        'species_specific': subset['human_specific'].sum() + subset['mouse_specific'].sum()
    })

pd.DataFrame(classification_summary).to_csv('risk_gene_classification_summary.csv', index=False)
