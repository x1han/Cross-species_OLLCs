"""
Single-Cell Data Integration using Pegasus and scVI

This module performs data integration across species using:
1. iNMF (Integrative Nonnegative Matrix Factorization) via Pegasus
2. scVI (single-cell Variational Inference) for OLLC subpopulation
"""

import pegasus as pg
import pegasusio as pgio
import scanpy as sc
import anndata as ad
import scvi
import numpy as np
import pandas as pd

# ============================================================================
# Full Atlas Integration using iNMF
# ============================================================================

# Load merged single-cell dataset
pgdata = pg.read_input('merged.h5ad')

# Identify robust genes for integration
pg.identify_robust_genes(pgdata)

# Apply log-normalization
pg.log_norm(pgdata)

# Perform iNMF integration with batch correction
inmf_key = pg.integrative_nmf(pgdata, batch = 'orig.ident', use_gpu = False)

# Build neighborhood graph using iNMF representation
pg.neighbors(pgdata, rep = inmf_key)

# Generate low-dimensional embeddings
pg.umap(pgdata, rep = inmf_key)
pg.tsne(pgdata, rep = inmf_key)

# Perform Louvain clustering
pg.louvain(pgdata, rep = inmf_key, resolution = 2)

# Save integrated atlas
pg.write_output(pgdata, 'atlas_full.h5ad')

# ============================================================================
# OLLC Subpopulation Analysis using scVI
# ============================================================================

# Subset to oligodendrocyte lineage cells (OPC and OL)
pgdata_ollcs = pgdata[pgdata.obs.cell_type.isin(['OL', 'OPC']), :].copy()

# Re-identify robust genes for subset
pg.identify_robust_genes(pgdata_ollcs)
pg.log_norm(pgdata_ollcs)

# Setup anndata for scVI model with batch and covariate information
pgdata_ollcs = scvi.model.SCVI.setup_anndata(pgdata_ollcs.to_anndata(), layer='X', batch_key='orig.ident', categorical_covariate_keys=['specie'])

# Train scVI model for latent representation
model = scvi.model.SCVI(pgdata_ollcs, n_hidden=128, n_latent=30, n_layers=2, gene_likelihood='nb')
model.train(max_epochs=400, early_stopping=True, early_stopping_patience=15, train_size=0.9)

# Extract latent representation from scVI
pgdata_ollcs.obsm['X_scvi'] = model.get_latent_representation()

# Build neighborhood graph on scVI embedding
sc.pp.neighbors(pgdata_ollcs, use_rep='X_scvi', n_neighbors=20)

# Generate UMAP visualization
sc.tl.umap(pgdata_ollcs)

# Perform Louvain clustering for cell state identification
sc.tl.louvain(pgdata_ollcs, resolution=1.5)

# Save OLLC-specific integrated dataset
pg.write_output(pgdata_ollcs, 'atlas_ollcs.h5ad')

