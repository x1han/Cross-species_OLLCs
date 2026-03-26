# ============================================================================
# Cross-Species Single-Cell Data Integration using SCTransform
# ============================================================================
#
# This script performs data integration across multiple species using:
# 1. SCTransform normalization for variance stabilization
# 2. Integration anchor identification for batch correction
# 3. Reference-based integration (Marques dataset as reference)
# ============================================================================

library(Seurat)

# ============================================================================
# SCTransform-based Integration Pipeline
# ============================================================================

# Load reference dataset for label transfer
marques_ref = trqwe::mcreadRDS('ollcs_label_transfer/ref/Marques.seurat.rds')

# Apply SCTransform normalization to all samples
for (i in names(seu.list)) {
    seu.list[[i]] <- SCTransform(seu.list[[i]], verbose = FALSE)
}

# Select integration features (3000 highly variable genes)
seu.features <- SelectIntegrationFeatures(object.list = seu.list, nfeatures = 3000)

# Prepare objects for SCT integration
seu.list <- PrepSCTIntegration(object.list = seu.list, anchor.features = seu.features)

# Designate reference dataset for integration
reference_dataset <- which(names(seu.list) == "SeuratProject")

# Find integration anchors between datasets
seu.anchors <- FindIntegrationAnchors(object.list = seu.list, normalization.method = "SCT", anchor.features = seu.features, reference = reference_dataset)

# Integrate data using identified anchors
seu.integrated <- IntegrateData(anchorset = seu.anchors, normalization.method = "SCT")

# Perform PCA dimensionality reduction
seu.integrated <- RunPCA(object = seu.integrated, verbose = FALSE)

# Generate UMAP embedding for visualization
seu.integrated <- RunUMAP(object = seu.integrated, dims = 1:20, n.neighbors = 50L)

# Visualize integrated data with cell type annotations
DimPlot(seu.integrated, group.by = 'cell_type', label = T)
