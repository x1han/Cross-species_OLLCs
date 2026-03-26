# ============================================================================
# Single-Cell RNA-seq Data Preprocessing Pipeline
# ============================================================================
# 
# This script provides preprocessing functions for single-cell data analysis:
# 1. Gene identifier conversion (ENSEMBL to official symbols)
# 2. Quality control and cell filtering
# 3. Mitochondrial gene percentage calculation
# 4. Doublet detection and removal
# 5. Cell cycle scoring
# ============================================================================

library(Seurat)
library(dplyr)

# ============================================================================
# Helper Functions for Data Preprocessing
# ============================================================================

#' Convert gene identifiers to official symbols using Ensembl database
#' @param taxid Species taxonomy ID
#' @return Processed transbase data frame with ENSEMBL and Synonym columns
make_syn_df = function(taxid) {
    # Read gene annotation database
    transbase <- read.delim(paste0('database/Enrichment/', taxid, '/', taxid, '.gene_info.txt'))
    # Extract ENSEMBL IDs from dbXrefs
    ENS <- unlist(lapply(strsplit(transbase$dbXrefs, split = '\\|'), function(x){x <- x[length(x)]}))
    ENS[setdiff(c(1:length(ENS)), grep(pattern = 'Ensembl', x = ENS))] <- '-'
    ENS_res <- unlist(lapply(strsplit(ENS, split = '\\:'), function(x){x <- x[length(x)]}))
    transbase$ENSEMBL <- ENS_res
    
    # Expand synonyms into separate rows
    transbase_res = data.frame(matrix('', length(unlist(strsplit(transbase$Synonyms, split = '[|]'))), ncol(transbase)))
    colnames(transbase_res) = colnames(transbase)
    j = 1
    for(i in 1:nrow(transbase)) {
        if(transbase$Synonyms[i] != '-') {
            synonyms_sp = strsplit(transbase$Synonyms[i], split = '[|]')[[1]]
            k = j+length(synonyms_sp)-1
            transbase_res[j:k, ] = transbase[i, ]
            transbase_res$Synonyms[j:k] = synonyms_sp
            j = k + 1
        } else {
            transbase_res[j, ] = transbase[i, ]
            j = j +1
        }
    }
    return(transbase_res)
}

#' Convert Ensembl IDs to gene symbols
#' @param genelist List of gene identifiers
#' @param transbase_res Reference database with symbol mappings
#' @return Gene list with converted symbols
ensembl_2_symbol = function(genelist, transbase_res) {
    result = genelist
    for(i in seq_along(genelist)) {
        g <- genelist[i]
        ens_matloc = which(transbase_res$ENSEMBL == g)
        symb_matloc = which(transbase_res$Symbol == g)
        if(length(ens_matloc) == 1 & length(symb_matloc) == 0) {
           result[i] = transbase_res$Symbol[ens_matloc]
        }
    }
    return(result)
}

#' Convert synonym names to official gene symbols
#' @param genelist List of gene names (may include synonyms)
#' @param transbase_res Reference database with symbol mappings
#' @return Gene list with converted symbols
synonyms_2_symbol = function(genelist, transbase_res) {
    result = genelist
    for(i in seq_along(genelist)) {
        g <- genelist[i]
        syno_matloc = which(transbase_res$Synonyms == g)
        symb_matloc = which(transbase_res$Symbol == g)
        if(length(syno_matloc) == 1 & length(symb_matloc) == 0) {
           result[i] = transbase_res$Symbol[syno_matloc]
        }
    }
    return(result)
}

#' Average expression values for duplicate genes
#' @param mtx Expression matrix with potentially duplicate gene names
#' @return Matrix with averaged duplicate genes
avg_dup_gene = function(mtx) {
    genelist = rownames(mtx)
    for(g in unique(genelist)) {
        matloc = which(rownames(mtx) == g)
        if(length(matloc) > 1) {
            mean_by_g = round(apply(mtx[matloc, ], 2, mean))
            mtx[matloc[1], ] = mean_by_g
            mtx = mtx[-matloc[2:length(matloc)], ]
        }
    }
    return(mtx)
}

#' Create one-to-one homologene mapping by removing duplicates
#' @param df Data frame with homologene pairs
#' @param check_col1 First column to check for duplicates
#' @param check_col2 Second column to check for duplicates
#' @return Filtered data frame with one-to-one mappings
homolo_make_one2one = function(df, check_col1, check_col2) {
    df = df[, c(check_col1, check_col2)]
    df = unique(df)
    # Remove duplicates from first column
    if(length(which(duplicated(df[, 1]))) != 0) {
        gene_dup = unique(df[which(duplicated(df[, 1])), 1])
        df = df[-which(df[, 1] %in% gene_dup), ]
    }
    # Remove duplicates from second column
    if(length(which(duplicated(df[, 2]))) != 0) {
        gene_dup = unique(df[which(duplicated(df[, 2])), 2])
        df = df[-which(df[, 2] %in% gene_dup), ]
    }
    return(data.frame(df))
}

#' Filter cells based on QC metrics
#' @param seu Seurat object
#' @return Filtered Seurat object
filter_cells = function(seu) {
    # Calculate log-transformed metrics
    seu$log10_nCount = log10(seu$nCount_RNA)
    seu$log10_nFeature = log10(seu$nFeature_RNA)
    
    # Calculate mean and standard deviation
    nCount_mean = mean(seu$log10_nCount)
    nCount_sd = sd(seu$log10_nCount)
    nFeature_mean = mean(seu$log10_nFeature)
    nFeature_sd = sd(seu$log10_nFeature)
    
    # Define thresholds (mean ± 1.96*SD)
    nCount_lower = 10^(nCount_mean - 1.96 * nCount_sd)
    nCount_upper = 10^(nCount_mean + 1.96 * nCount_sd)
    nFeature_lower = 10^(nFeature_mean - 1.96 * nFeature_sd)
    nFeature_upper = 10^(nFeature_mean + 1.96 * nFeature_sd)
    
    # Apply filters: UMIs, genes, and mitochondrial percentage
    keep_cells = seu$nCount_RNA >= nCount_lower & 
                 seu$nCount_RNA <= nCount_upper &
                 seu$nFeature_RNA >= nFeature_lower & 
                 seu$nFeature_RNA <= nFeature_upper &
                 seu$percent_mito < 0.2
    
    seu = subset(seu, cells = Cells(seu)[keep_cells])
    return(seu)
}

#' Process single-cell data for a specific species
#' @param seu_list List of Seurat objects
#' @param taxid Species taxonomy ID
#' @return Merged and processed Seurat object
process_species = function(seu_list, taxid) {
    # Split large samples (>10k cells) for efficient processing
    seu_list_split = lapply(seu_list, function(seu) {
        ncell = ncol(seu)
        if (ncell > 10000) {
            n = ceiling(ncell/10000)
            seu$id = sample(1:n, ncell, replace = T)
            seu = SplitObject(seu, "id")
        }
        return(seu)
    })
    
    seu_list_split = unlist(seu_list_split)
    transbase_res = make_syn_df(taxid)
    
    # Convert gene identifiers and filter cells
    seu_list_split = lapply(seu_list_split, function(seu) {
        counts_raw = seu@assays$RNA@counts
        meta_raw = seu@meta.data
        features_raw = rownames(seu)
        
        # Convert Ensembl and synonyms to official symbols
        features_new = ensembl_2_symbol(features_raw, transbase_res = transbase_res)
        features_new = synonyms_2_symbol(features_new, transbase_res = transbase_res)
        
        counts_new = counts_raw
        rownames(counts_new) = features_new
        counts_new = avg_dup_gene(counts_new)
        
        seu_new = CreateSeuratObject(counts_new)
        seu_new@meta.data = meta_raw
        seu_new = filter_cells(seu_new)
        
        return(seu_new)
    })
    
    # Find intersection genes across all samples
    feature_list = lapply(seu_list_split, rownames)
    inter_gene = Reduce(intersect, feature_list)
    
    # Subset to common genes and merge
    seu_list = lapply(seu_list_split, function(seu) {
        subset(seu, features = inter_gene)
    })
    
    seu_merged = Reduce(merge, seu_list)
    return(seu_merged)
}

# ============================================================================
# Homologene Mapping Construction
# ============================================================================

# Load biomaRt databases for cross-species mapping
mart <- biomaRt::useEnsembl(biomart = "ensembl", mirror = "asia")
human = readRDS('biomaRt/human.rds')
mouse = readRDS('biomaRt/mouse.rds')
macaca_fasci = readRDS('biomaRt/macaca_fasci.rds')
macaca_mulatta = readRDS('biomaRt/macaca_mulatta.rds')

# Generate mouse-human ortholog mapping
mouse.gene = rownames(mouse_seu)
mouse2h.g <- biomaRt::getLDS(attributes = c("entrezgene_id", "ensembl_gene_id", "mgi_symbol"), filters = "mgi_symbol",
                values = mouse.gene, mart = mouse,
                attributesL = c("entrezgene_id", "ensembl_gene_id","hgnc_symbol"),
                martL = human, uniqueRows = T)

# Create one-to-one mappings for all species pairs
homolo_one2one_mouse = homolo_make_one2one(mouse2h.g, 'MGI.symbol', 'HGNC.symbol')
homolo_one2one_macaca_fasci = homolo_make_one2one(macaca_fasci2h.g, 'HGNC.symbol', 'HGNC.symbol.1')
homolo_one2one_macaca_mulatta = homolo_make_one2one(macaca_mulatta2h.g, 'HGNC.symbol', 'HGNC.symbol.1')

# Merge all homologene mappings
homolo_list = list(homolo_one2one_mouse, homolo_one2one_macaca_fasci, homolo_one2one_macaca_mulatta)
homolo_df = homolo_list[[1]]
for(i in 2:length(homolo_list)) {
    homolo_df = merge(homolo_df, homolo_list[[i]], by = colnames(homolo_list[[i]])[1])
}
colnames(homolo_df) = c('human', 'mouse', 'macaca_fascicularis', 'macaca_mulatta')

# Ensure strict one-to-one relationships across all four species
homolo_df = homolo_make_one2one(homolo_df, 'human', 'mouse')
homolo_df = homolo_make_one2one(homolo_df, 'human', 'macaca_fascicularis')
homolo_df = homolo_make_one2one(homolo_df, 'human', 'macaca_mulatta')

# Save homologene table
write.table(homolo_df, 'homologene/homologene.xls', sep = '\t', quote = F, row.names = F)

# ============================================================================
# Alternative Homologene Sources Comparison
# ============================================================================

# Load homologene mappings from different sources
mouse2human.biomart = read.delim('biomaRt/mouse2human.homologene.xls')
mouse2human.naka = read.delim('Nakamura_et_al/mouse2human.homologene.xls')
fasci2human.biomart = read.delim('biomaRt/macaca_fasci2human.homologene.xls')
fasci2human.naka = read.delim('Nakamura_et_al/macaca_fasci2human.homologene.xls')
mulatta2human.biomart = read.delim('biomaRt/macaca_mulatta2human.homologene.xls')

# Convert all to one-to-one mappings for comparison
mouse2human.biomart = homolo_make_one2one(mouse2human.biomart, 2, 4)
mouse2human.naka = homolo_make_one2one(mouse2human.naka, 2, 4)
fasci2human.biomart = homolo_make_one2one(fasci2human.biomart, 2, 4)
fasci2human.naka = homolo_make_one2one(fasci2human.naka, 2, 4)
mulatta2human.biomart = homolo_make_one2one(mulatta2human.biomart, 2, 4)

# ============================================================================
# Load Single-Cell Data for Processing
# ============================================================================

# Load final homologene mapping
homolo = read.delim('homologene_human_mouse_fasci_mulatta.xls')

# Load Seurat objects for all species
human_seu = trqwe::mcreadRDS('human.seurat.rds', mc.cores = 32)
mouse_seu = trqwe::mcreadRDS('mouse.seurat.rds', mc.cores = 32)
fasci_seu = trqwe::mcreadRDS('macaca_fascicularis.seurat.rds', mc.cores = 32)
mulatta_seu = trqwe::mcreadRDS('macaca_mulatta.seurat.rds', mc.cores = 32)

# ============================================================================
# Quality Control Visualization
# ============================================================================

# Plot QC metrics distribution across samples
QC_df = rbind(
    data.frame(nCount_RNA = human_seu$nCount_RNA, nFeature_RNA = human_seu$nFeature_RNA, source = 'human'),
    data.frame(nCount_RNA = mouse_seu$nCount_RNA, nFeature_RNA = mouse_seu$nFeature_RNA, source = 'mouse')
)

p1 = ggplot(QC_df, aes(x = log10(nCount_RNA), color = source, fill = source)) + 
     geom_histogram(binwidth = 0.01) + 
     geom_density(alpha = 0.5) + 
     scale_fill_manual(values = collect_cols(1:2)) + 
     scale_color_manual(values = collect_cols(1:2)) + 
     ggtitle('log10(nCount_RNA)') + theme_bw()
p2 = ggplot(QC_df, aes(x = log10(nFeature_RNA), color = source, fill = source)) + 
     geom_histogram(binwidth = 0.01) + 
     geom_density(alpha = 0.5) + 
     scale_fill_manual(values = collect_cols(1:2)) + 
     scale_color_manual(values = collect_cols(1:2)) + 
     ggtitle('log10(nFeature_RNA)') + theme_bw()

p = p1/p2


