# ============================================================================
# Pseudobulk Differential Expression Analysis Pipeline
# ============================================================================
#
# This script performs pseudobulk-based differential expression analysis using:
# 1. DESeq2 for statistical modeling and hypothesis testing
# 2. RUVr (Remove Unwanted Variation) for batch effect correction
# 3. LFC shrinkage for improved effect size estimation
#
# Key Features:
# - Aggregates single-cell counts to sample-level pseudobulks
# - Handles complex experimental designs with covariates
# - Generates comprehensive QC plots and result visualizations
# - Supports multiple species and cell types
# ============================================================================

suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(RUVSeq))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(FactoMineR))
suppressPackageStartupMessages(library(limma))

set.seed(492299)

# ============================================================================
# CONFIGURATION - Modify these parameters before running
# ============================================================================

# Input/Output paths
SEURAT_RDS <- 'seurat_object.rds'           # Path to input Seurat object
OUTPUT_DIR <- 'DESeq2_results'              # Output directory

# Column names in meta.data
BY_COL <- 'group'                           # Grouping variable column name
SAMPLE_COL <- 'sample'                      # Sample ID column for pseudobulk aggregation

# Comparison groups
CASE_NAME <- 'case'                         # Case/treatment group name
REF_NAME <- 'ref'                           # Reference/control group name

# Analysis options
DO_PSEUDOBULK <- TRUE                       # If TRUE, aggregate cells to pseudobulk
REMOVE_BATCH <- TRUE                        # If TRUE, perform RUVr batch correction
MIN_EXP <- 10                               # Minimum expression threshold
PVALUE_THRESHOLD <- 0.05                    # P-value significance threshold
LOG2FC_THRESHOLD <- 1                       # Log2FC threshold
TAXID <- 9606                               # Taxonomy ID (9606=human, 10090=mouse, etc.)
USE_SHRINK <- TRUE                          # If TRUE, perform LFC shrinkage
USE_LRT <- FALSE                            # If TRUE, use LRT test instead of Wald

# Cell type/subtype subsetting (optional)
SUBSET_BY <- NULL                           # Column name for subsetting (e.g., 'cell_type')
SUBSET_VALUE <- NULL                        # Value(s) to keep (comma-separated if multiple)

# ============================================================================
# Helper Functions
# ============================================================================

#' Convert Ensembl IDs to gene symbols
#' @param genelist Vector of gene identifiers
#' @param transbase_res Data frame with gene annotation
#' @return Vector of gene symbols
ensembl_2_symbol <- function(genelist, transbase_res) {
    result <- genelist
    for(i in seq_along(genelist)) {
        g <- genelist[i]
        ens_matloc <- which(transbase_res$ENSEMBL == g)
        symb_matloc <- which(transbase_res$Symbol == g)
        if(length(ens_matloc) == 1 & length(symb_matloc) == 0) {
           result[i] <- transbase_res$Symbol[ens_matloc]
        }
    }
    return(result)
}

#' Convert synonyms to gene symbols
#' @param genelist Vector of gene identifiers
#' @param transbase_res Data frame with gene annotation
#' @return Vector of gene symbols
synonyms_2_symbol <- function(genelist, transbase_res) {
    result <- genelist
    for(i in seq_along(genelist)) {
        g <- genelist[i]
        syno_matloc <- which(transbase_res$Synonyms == g)
        symb_matloc <- which(transbase_res$Symbol == g)
        if(length(syno_matloc) == 1 & length(symb_matloc) == 0) {
           result[i] <- transbase_res$Symbol[syno_matloc]
        }
    }
    return(result)
}

#' Average duplicate genes
#' @param mtx Expression matrix
#' @return Matrix with averaged duplicates
avg_dup_gene <- function(mtx) {
    genelist <- rownames(mtx)
    for(g in unique(genelist)) {
        matloc <- which(rownames(mtx) == g)
        if(length(matloc) > 1) {
            mean_by_g <- round(apply(mtx[matloc, ], 2, mean))
            mtx[matloc[1], ] <- mean_by_g
            mtx <- mtx[-matloc[2:length(matloc)], ]
        }
    }
    return(mtx)
}

# ============================================================================
# MAIN PIPELINE
# ============================================================================

# Create output directory
if(!dir.exists(OUTPUT_DIR)){
    dir.create(OUTPUT_DIR, recursive = TRUE)
}

# --- 1. Data Loading and Preprocessing ---
seu <- readRDS(SEURAT_RDS)

# Subset cells if specified
if(!is.null(SUBSET_BY)){
    subvalue_vec <- strsplit(SUBSET_VALUE, split = ',')[[1]]
    seu <- subset(seu, cells = Cells(seu)[which(seu@meta.data[, SUBSET_BY] %in% subvalue_vec)])
}

# Keep only case and ref groups
seu <- subset(seu, cells = Cells(seu)[which(seu@meta.data[, BY_COL] %in% c(CASE_NAME, REF_NAME))])

metadata <- seu@meta.data
counts <- seu@assays$RNA@counts

# --- 2. Pseudobulk Construction ---
if(DO_PSEUDOBULK) {
    samp_fac <- factor(metadata[, SAMPLE_COL])
    mm <- model.matrix(~ 0 + samp_fac)
    colnames(mm) <- levels(samp_fac)
    countData <- counts %*% mm
    countData <- as.matrix(countData)
} else {
    countData <- as.matrix(counts)
}

# Prepare sample metadata
unique_meta <- unique(metadata[, c(SAMPLE_COL, BY_COL)])
rownames(unique_meta) <- unique_meta[, SAMPLE_COL]
colData <- unique_meta[colnames(countData), ]
colData <- data.frame(
  row.names = rownames(colData), 
  group = factor(colData[, BY_COL], levels = c(REF_NAME, CASE_NAME))
)

# --- 3. Filter Lowly Expressed Genes ---
keep_rows <- apply(countData, 1, function(row) {
    sum(row > MIN_EXP) >= (length(row) * 0.5)
})
countData <- countData[keep_rows, , drop = FALSE]

# --- 4. Batch Effect Estimation (RUVr with LVGs) ---
design_formula <- ~ group
reduced_formula <- ~ 1 

if(REMOVE_BATCH) {
    # Step 4.1: Identifying Lowly Variable Genes (LVGs) as empirical controls...
    
    # Create temporary object to calculate standardized variance
    tmp_seu <- CreateSeuratObject(counts = countData)
    tmp_seu <- NormalizeData(tmp_seu, verbose = FALSE)
    tmp_seu <- FindVariableFeatures(tmp_seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    
    # Extract genes with lowest standardized variance as control genes
    lvg_info <- tmp_seu[["RNA"]]@meta.features %>%
        as.data.frame() %>%
        filter(vst.mean > 0.1) %>%
        arrange(vst.variance.standardized)
    
    lvg_list <- head(rownames(lvg_info), 2000)
    lvg_list <- intersect(lvg_list, rownames(countData))

    seq <- newSeqExpressionSet(as.matrix(round(countData)))
    design_temp <- model.matrix(~group, data=colData)
    dge <- DGEList(counts = counts(seq), group = colData$group)
    dge <- calcNormFactors(dge, method="RLE")
    dge <- estimateGLMCommonDisp(dge, design_temp)
    dge <- estimateGLMTagwiseDisp(dge, design_temp)
    fit <- glmFit(dge, design_temp)
    res <- residuals(fit, type="deviance")
    
    seqUQ <- betweenLaneNormalization(seq, which="upper")
    
    # Use lvg_list as control genes (cIdx)
    seqRUVr <- RUVr(seqUQ, lvg_list, k=1, res)
    
    # Extract W1 and add to colData
    colData$W1 <- pData(seqRUVr)$W_1
    
    design_formula <- ~ W1 + group
    reduced_formula <- ~ W1 
    
} else {
    # Skip batch correction
}

# --- 5. Run DESeq2 ---
dds <- DESeqDataSetFromMatrix(countData = round(countData), colData = colData, design = design_formula)
dds$group <- relevel(dds$group, ref = REF_NAME)

if(USE_LRT) {
    dds_norm <- DESeq(dds, test = "LRT", reduced = reduced_formula)
} else {
    dds_norm <- DESeq(dds)
}

# --- 6. Quality Control Plots ---

# Perform VST transformation
vsd <- varianceStabilizingTransformation(dds_norm, blind = FALSE)
vsd_mat <- assay(vsd)

# Remove batch effect from visualization if batch correction was performed
if(REMOVE_BATCH) {
    vsd_mat <- limma::removeBatchEffect(vsd_mat, covariates = colData(dds_norm)[["W1"]])
}

# 6.1 PCA
pca_res <- PCA(t(vsd_mat), graph = FALSE)
pca_plot <- fviz_pca_ind(pca_res,
                        geom = c("point"), 
                        col.ind = colData$group, 
                        palette = c("#00AFBB", "#E7B800"),
                        addEllipses = TRUE,
                        legend.title = "Group",
                        title = ifelse(REMOVE_BATCH, "PCA (Batch Removed)", "PCA (Raw)")) + theme_bw()

ggsave(plot = pca_plot, paste0(OUTPUT_DIR, '/PCA.png'), dpi = 300)
ggsave(plot = pca_plot, paste0(OUTPUT_DIR, '/PCA.pdf'))
pdf(paste0(OUTPUT_DIR, '/PCA_dims.pdf'))
factoextra::fviz_eig(pca_res, addlabels = T, ylim = c(0, 100))
dev.off()

# 6.2 Hierarchical Clustering
sampleDists <- dist(t(vsd_mat))
hc <- hclust(sampleDists, method = "ward.D2")
pdf(paste0(OUTPUT_DIR, '/hclust.pdf'))
plot(hc, hang = -1, main = ifelse(REMOVE_BATCH, "Cluster Dendrogram (Batch Removed)", "Cluster Dendrogram"))
dev.off()

# --- 7. Extract Results with Shrinkage ---
contrast_vec <- c('group', CASE_NAME, REF_NAME)
res_raw <- results(dds_norm, contrast = contrast_vec, alpha = PVALUE_THRESHOLD)

pdf(paste0(OUTPUT_DIR, '/MAplot_without_shrink.pdf'))
plotMA(res_raw, ylim=c(-4,4), main="Without Shrinkage")
dev.off()

if(USE_SHRINK){
    # Perform LFC Shrinkage
    # Fallback logic for missing ashr
    if (requireNamespace("ashr", quietly = TRUE)) {
        res_shrink <- lfcShrink(dds_norm, contrast = contrast_vec, res = res_raw, type = "ashr")
    } else {
        res_shrink <- lfcShrink(dds_norm, contrast = contrast_vec, res = res_raw, type = "normal")
    }
    
    pdf(paste0(OUTPUT_DIR, '/MAplot_after_shrink.pdf'))
    plotMA(res_shrink, ylim=c(-4,4), main="After Shrinkage")
    dev.off()
    dd_res <- res_shrink
} else {
    dd_res <- res_raw
}

# --- 8. Process and Export Results ---
res_df_raw <- as.data.frame(dd_res) %>% 
    na.omit() %>% 
    dplyr::filter(!is.na(log2FoldChange)) %>% 
    dplyr::mutate(FoldChange = 2^log2FoldChange) %>% 
    dplyr::arrange(desc(log2FoldChange))

res_df_raw$up_down <- ifelse(res_df_raw$log2FoldChange > 0, "Up", "Down")

# Note: Gene annotation conversion requires transbase function (if available)
# Uncomment if you have the transbase function
# if(exists("transbase")){
#     res_df_raw$gene <- tryCatch({
#         transbase(genelist = rownames(res_df_raw), taxid = TAXID, from = 'ENSEMBL', to = 'Symbol')
#     }, error = function(e) { rownames(res_df_raw) }) 
#     
#     res_df_raw <- cbind(res_df_raw, data.frame(counts(dds_norm, normalized=TRUE)[rownames(res_df_raw), ], check.names=F))
#     
#     res_df_raw <- tryCatch({
#         geneanno(res_df_raw, taxid = TAXID)
#     }, error = function(e) { res_df_raw })
# }

# For now, use rownames as gene names
res_df_raw$gene <- rownames(res_df_raw)
res_df_raw <- cbind(res_df_raw, data.frame(counts(dds_norm, normalized=TRUE)[rownames(res_df_raw), ], check.names=F))

# Filter significant results
res_sig_p <- res_df_raw %>% dplyr::filter(pvalue < PVALUE_THRESHOLD)
res_sig_fc <- res_df_raw %>% dplyr::filter(pvalue < PVALUE_THRESHOLD, abs(log2FoldChange) > LOG2FC_THRESHOLD)

# Define columns to save
cols_to_save <- c("gene", "pvalue", "padj", "FoldChange", "log2FoldChange", "baseMean", "up_down")
if("counts" %in% colnames(res_df_raw)){
    count_cols <- grep("^\\d+", colnames(res_df_raw), value = TRUE)
    cols_to_save <- c(cols_to_save, count_cols)
}
ideal_order <- unique(c(cols_to_save, "description", "GO_term"))
existing_cols <- intersect(ideal_order, colnames(res_df_raw))

# Export results
write.table(res_df_raw[, existing_cols], 
            paste0(OUTPUT_DIR, '/DESeq2.result.raw.xls'), 
            sep = '\t', quote = F, row.names = F)
write.table(res_sig_p[, existing_cols], 
            paste0(OUTPUT_DIR, '/DESeq2.result.filt_pvalue.xls'), 
            sep = '\t', quote = F, row.names = F)
write.table(res_sig_fc[, existing_cols], 
            paste0(OUTPUT_DIR, '/DESeq2.result.filtered.xls'), 
            sep = '\t', quote = F, row.names = F)

# --- 9. Heatmap Visualization ---
top_up <- res_sig_fc %>% filter(up_down == 'Up') %>% top_n(15, log2FoldChange) %>% pull(gene)
top_down <- res_sig_fc %>% filter(up_down == 'Down') %>% top_n(15, -log2FoldChange) %>% pull(gene)
top_genes <- c(top_up, top_down)

target_ids <- rownames(res_df_raw)[res_df_raw$gene %in% top_genes]
target_ids <- target_ids[target_ids %in% rownames(vsd_mat)]

if(length(target_ids) > 2) {
    mat_heatmap <- vsd_mat[target_ids, ]
    match_idx <- match(rownames(mat_heatmap), rownames(res_df_raw))
    rownames(mat_heatmap) <- res_df_raw$gene[match_idx]
    
    anno_col <- data.frame(group = colData(dds_norm)$group)
    rownames(anno_col) <- colnames(dds_norm)
    
    p <- pheatmap(mat_heatmap, cluster_rows = TRUE, cluster_cols = TRUE, 
                  show_rownames = TRUE, annotation_col = anno_col, 
                  scale = "row", main = paste0("Top DE Genes: ", CASE_NAME, " vs ", REF_NAME))
    
    ggsave(plot = p, paste0(OUTPUT_DIR, '/top_Diffexp_heatmap.pdf'))
    ggsave(plot = p, paste0(OUTPUT_DIR, '/top_Diffexp_heatmap.png'), dpi = 300)
} else {
    # Skip heatmap if insufficient DEGs
}

# --- 10. Volcano Plot ---
df_vo <- res_df_raw
df_vo$Regulated <- 'ns'
df_vo$Regulated[df_vo$log2FoldChange >= LOG2FC_THRESHOLD & df_vo$pvalue < PVALUE_THRESHOLD] <- 'up'
df_vo$Regulated[df_vo$log2FoldChange <= -LOG2FC_THRESHOLD & df_vo$pvalue < PVALUE_THRESHOLD] <- 'down'

top_n_label <- 10
label_data <- df_vo %>%
    dplyr::filter(Regulated != 'ns') %>%
    dplyr::group_by(Regulated) %>%
    dplyr::arrange(pvalue) %>% 
    dplyr::slice_head(n = top_n_label) %>%
    dplyr::ungroup()

p <- ggplot(data = df_vo, aes(x = log2FoldChange, y = -log10(pvalue), color = Regulated)) +
    geom_point(size = 1) +
    scale_color_manual(values = c('red', 'gray', 'blue'), limits = c('up', 'ns', 'down')) +
    labs(x = 'log2 FoldChange', y = '-log10 p-value') +
    geom_vline(xintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD), lty = 3, color = 'black') +
    geom_hline(yintercept = -log10(PVALUE_THRESHOLD), lty = 3, color = 'black') + 
    geom_text_repel(data = label_data, aes(label = gene), 
                    arrow = arrow(length=unit(0.01, "npc")), 
                    show.legend = F, 
                    box.padding = 0.5,
                    max.overlaps = Inf) + 
    theme_classic() + 
    guides(fill = guide_legend(title = '', override.aes = aes(label = '')))

ggsave(plot = p, paste0(OUTPUT_DIR, '/Diffexp_', CASE_NAME, '-v-', REF_NAME, '_volcano.pdf'))
ggsave(plot = p, paste0(OUTPUT_DIR, '/Diffexp_', CASE_NAME, '-v-', REF_NAME, '_volcano.png'), dpi = 300)

# ============================================================================
# COMPLETION
# ============================================================================
