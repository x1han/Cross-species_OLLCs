"""
Cross-Species Cell Type Classification using Machine Learning

This module evaluates species-specific cell type signatures using:
1. MLP (Multi-Layer Perceptron)
2. XGBoost (Extreme Gradient Boosting)
3. LightGBM (Light Gradient Boosting Machine)
"""

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import accuracy_score, confusion_matrix
import xgboost as xgb
import lightgbm as lgb
import optuna
import matplotlib.pyplot as plt
import seaborn as sns

# ============================================================================
# Load Data and Prepare Features
# ============================================================================

# Load integrated single-cell dataset
adata = sc.read_h5ad('integrated.h5ad')

# Load pseudobulk differential expression results
de_results = pd.read_csv('pseudobulk_DE_results.csv')

# Extract significant DEGs (adjusted p-value < 0.05)
de_genes = de_results[de_results['padj'] < 0.05]['gene'].tolist()

# Extract expression matrix for DEGs only
X = adata[:, de_genes].X
if hasattr(X, 'toarray'):
    X = X.toarray()

# Get cell type labels
y = adata.obs['cell_type'].values

# Normalize features to [0, 1] range
scaler = MinMaxScaler()
X_scaled = scaler.fit_transform(X)

# ============================================================================
# Cross-Species Classification Pipeline
# ============================================================================

# Get unique species for leave-one-species-out validation
species_list = adata.obs['specie'].unique()
results = []

# Iterate through each species as training set
for train_species in species_list:
    # Create train/test split by species
    train_mask = adata.obs['specie'] == train_species
    test_mask = ~train_mask
    
    X_train = X_scaled[train_mask]
    y_train = y[train_mask]
    X_test = X_scaled[test_mask]
    y_test = y[test_mask]
    
    # Train MLP classifier with optimized hyperparameters
    mlp = MLPClassifier(hidden_layer_sizes=(180,), activation='relu', 
                       batch_size='auto', learning_rate_init=1e-5,
                       max_iter=500, early_stopping=True, 
                       validation_fraction=0.2, n_iter_no_change=15,
                       random_state=492299)
    mlp.fit(X_train, y_train)
    y_pred_mlp = mlp.predict(X_test)
    acc_mlp = accuracy_score(y_test, y_pred_mlp)
    
    # Train XGBoost classifier
    xgb_model = xgb.XGBClassifier(n_estimators=100, max_depth=6, 
                                  learning_rate=0.1, random_state=492299,
                                  use_label_encoder=False, eval_metric='mlogloss')
    xgb_model.fit(X_train, y_train)
    y_pred_xgb = xgb_model.predict(X_test)
    acc_xgb = accuracy_score(y_test, y_pred_xgb)
    
    # Train LightGBM classifier
    lgb_model = lgb.LGBMClassifier(n_estimators=100, max_depth=6,
                                   learning_rate=0.1, random_state=492299)
    lgb_model.fit(X_train, y_train)
    y_pred_lgb = lgb_model.predict(X_test)
    acc_lgb = accuracy_score(y_test, y_pred_lgb)
    
    # Store classification results
    results.append({
        'train_species': train_species,
        'mlp_accuracy': acc_mlp,
        'xgboost_accuracy': acc_xgb,
        'lightgbm_accuracy': acc_lgb
    })

# Save cross-species classification performance metrics
results_df = pd.DataFrame(results)
results_df.to_csv('cross_species_classification.csv', index=False)
