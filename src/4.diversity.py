import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# ============================================================================
# OLLC Subtype Diversity Analysis (Based on Pseudobulk DESeq2 Results)
# ============================================================================

adata = sc.read_h5ad('integrated.h5ad')

# ============================================================================
# 1. Read Pseudobulk DESeq2 Results (Cross-species Comparison)
# ============================================================================

de_results = {}
stages = ['OPC', 'COP', 'NFO', 'MOL']
species_list = ['human', 'fasci', 'mulatta', 'mouse']

for stage in stages:
    de_results[stage] = {}
    for sp in species_list:
        # Read from pseudobulk result files
        file_path = f'specie_divergence/stage/pseudobulks/{stage}/group1/{sp}_vs_human/DESeq2.result.filtered.xls'
        
        if os.path.exists(file_path):
            de_df = pd.read_excel(file_path)
            
            # Extract significantly differentially expressed genes (|log2FC| > 1, adj.P < 0.05)
            deg_up = de_df[de_df['log2FoldChange'] > 0]['gene'].tolist()
            deg_down = de_df[de_df['log2FoldChange'] < 0]['gene'].tolist()
            
            de_results[stage][sp] = {
                'deg_up': set(deg_up),
                'deg_down': set(deg_down),
                'all_deg': set(de_df['gene'].tolist()),
                'full_results': de_df
            }
        else:
            de_results[stage][sp] = {'deg_up': set(), 'deg_down': set(), 'all_deg': set(), 'full_results': None}

# ============================================================================
# 2. Circos Plot Data Preparation (Conserved and Species-specific DEGs)
# ============================================================================

# Calculate conserved and species-specific genes for each stage
circos_data = []
for stage in stages:
    # Get DEGs from all species
    all_degs = set()
    for sp in species_list:
        if de_results[stage][sp]['all_deg']:
            all_degs.update(de_results[stage][sp]['all_deg'])
    
    # For each gene, mark in which species it is differentially expressed
    for gene in all_degs:
        present_in = [sp for sp in species_list if gene in de_results[stage][sp]['all_deg']]
        
        circos_data.append({
            'stage': stage,
            'gene': gene,
            'present_in_species': present_in,
            'num_species': len(present_in),
            'conserved': len(present_in) == len(species_list),  # Present in all species
            'primate_specific': all(sp in present_in for sp in ['human', 'fasci', 'mulatta']) and 'mouse' not in present_in,
            'mouse_specific': present_in == ['mouse'],
            'species_specific': len(present_in) == 1
        })

circos_df = pd.DataFrame(circos_data)
circos_df.to_csv('ollc_stage_degs_circos.csv', index=False)

# ============================================================================
# 3. PAGA Trajectory Inference (Keep Original Logic)
# ============================================================================

for stage in stages:
    adata_stage = adata[adata.obs['cell_type_stage'] == stage, :].copy()
    
    if len(adata_stage) == 0:
        continue
    
    sc.pp.neighbors(adata_stage, use_rep='X_scvi', n_neighbors=20)
    sc.tl.paga(adata_stage, groups='specie')
    
    # Visualize PAGA graph
    sc.pl.paga(adata_stage, threshold=0.1, node_size_scale=1.5, edge_width_scale=0.5, show=False)
    plt.savefig(f'paga_{stage}.pdf', bbox_inches='tight')
    plt.savefig(f'paga_{stage}.png', dpi=600, bbox_inches='tight')
    plt.close()

# ============================================================================
# 4. SCENIC Analysis Preparation (Keep Original Logic with Annotations)
# ============================================================================

# Prepare expression matrix for each species
for spec in adata.obs['specie'].unique():
    adata_spec = adata[adata.obs['specie'] == spec, :].copy()
    
    if len(adata_spec) == 0:
        continue
    
    # Convert expression matrix to loom format (for SCENIC)
    expr_mat = pd.DataFrame(
        adata_spec.X.toarray() if hasattr(adata_spec.X, 'toarray') else adata_spec.X,
        index=adata_spec.obs_names,
        columns=adata_spec.var_names
    )
    
    # Save as loom format (requires loompy)
    try:
        import loompy as lp
        lp.create(
            f'{spec}_expression.loom',
            {'Gene': expr_mat.values.T},
            {'CellID': expr_mat.index.values},
            {'RowAttrs': ['Gene']}
        )
    except ImportError:
        # Save as CSV as backup
        expr_mat.to_csv(f'{spec}_expression.csv')

# ============================================================================
# 5. Statistics for Conserved and Species-specific DEG Proportions
# ============================================================================

stats_summary = []
for stage in stages:
    total_degs = len(circos_df[circos_df['stage'] == stage])
    conserved_count = circos_df[(circos_df['stage'] == stage) & (circos_df['conserved'])].shape[0]
    primate_specific_count = circos_df[(circos_df['stage'] == stage) & (circos_df['primate_specific'])].shape[0]
    mouse_specific_count = circos_df[(circos_df['stage'] == stage) & (circos_df['mouse_specific'])].shape[0]
    species_specific_count = circos_df[(circos_df['stage'] == stage) & (circos_df['species_specific'])].shape[0]
    
    stats_summary.append({
        'stage': stage,
        'total_degs': total_degs,
        'conserved_count': conserved_count,
        'conserved_pct': conserved_count / total_degs * 100 if total_degs > 0 else 0,
        'primate_specific_count': primate_specific_count,
        'primate_specific_pct': primate_specific_count / total_degs * 100 if total_degs > 0 else 0,
        'mouse_specific_count': mouse_specific_count,
        'mouse_specific_pct': mouse_specific_count / total_degs * 100 if total_degs > 0 else 0,
        'species_specific_count': species_specific_count,
        'species_specific_pct': species_specific_count / total_degs * 100 if total_degs > 0 else 0
    })

stats_df = pd.DataFrame(stats_summary)
stats_df.to_csv('ollc_stage_deg_statistics.csv', index=False)
