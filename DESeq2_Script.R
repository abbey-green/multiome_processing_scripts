#!/usr/bin/env Rscript

# Load required libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(apeglm)
library(tidyverse)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
    stop("Usage: Rscript DESeq2_Script.R <region> <cell_type_list> <input_dir> <pseudobulk_dir> <output_dir>")
}

region <- args[1]
input_dir <- args[2]
pb_date <- args[3]
output_dir <- args[4]

# # Capture all remaining arguments as cell types
# cell_type_list <- args[5:length(args)]

# Print input parameters
print(paste("Region:", region))
print(paste("Input directory:", input_dir))
print(paste("Output directory:", output_dir))
# print(paste("Cell types:", paste(cell_type_list, collapse = ", ")))

PFC_cell_list <- list('Oligodendrocyte','L5 IT','L2-3 IT', 'Vip',
 'Pvalb', 'L6 IT Car3', 'L4 IT', 'Lamp5', 'Astrocyte', 'Lamp5 Lhx6', 'Microglia-PVM',
 'OPC','L5-6 NP','Sncg','Sst', 'L6b','L6 IT Car3', 'L6 IT', 'VLMC', 'L6 CT', 'Pax6','Endothelial', 'Chandelier')
AMYG_cell_list <- cell_list <- list('Oligodendrocyte', 'Astrocyte','Microglia','OPC','Endothelial','Fibroblast')
CB_cell_list <- list('Granule','Bergmann','OGD','OGP','Microglia', 'MLI Type I', 'Astrocyte','MLI Type II', 'Fibroblast')

conditions <- c('HIV', 'OUD', 'HIVOUD_Interaction')

if (region == 'PFC') {
    cell_type_list <- PFC_cell_list
    }

if (region == 'AMYG') {
    cell_type_list <- AMYG_cell_list
    }

if (region == 'CB') {
    cell_type_list <- CB_cell_list
    }

for (condition in conditions) {
    for (cell_type in cell_type_list) {
        output_dir_2 <- file.path(output_dir, condition)
        dir.create(output_dir_2, showWarnings = FALSE, recursive = TRUE)
        
        col_data <- read.csv(paste(input_dir,cell_type,'_metadata_',region,'_',pb_date,'.csv', sep = ""), row.names=1)
        count_data <- read.csv(paste(input_dir,cell_type,'_psuedobulk_',region,'_',pb_date,'.csv', sep = ""), 
                               header=TRUE, check.names=FALSE, row.names=1)
        
        # Ensure column names in count matrix match row names in metadata
        stopifnot(identical(colnames(count_data), rownames(col_data)))
    
        # Remove outlier CC108
        col_data <- col_data[rownames(col_data) != "CC108", ]
        count_data <- count_data[, !(colnames(count_data) == "CC108")]
    
        if (condition == 'HIV') {
            # Remove OUD and HIV+OUD donors
            hivoud_donors <- rownames(col_data)[col_data$condition == "HIV+OUD"]
            oud_donors <- rownames(col_data)[col_data$condition == "OUD"]
            col_data <- col_data[!rownames(col_data) %in% c(hivoud_donors, oud_donors), ]
            count_data <- count_data[, !colnames(count_data) %in% c(hivoud_donors, oud_donors)]
        }
    
        if (condition == 'OUD') {
            hivoud_donors <- rownames(col_data)[col_data$condition == "HIV+OUD"]
            hiv_donors <- rownames(col_data)[col_data$condition == "HIV"]
            col_data <- col_data[!rownames(col_data) %in% c(hivoud_donors, hiv_donors), ]
            count_data <- count_data[, !colnames(count_data) %in% c(hivoud_donors, hiv_donors)]
        }
    
        if (condition == 'HIVOUD_Interaction') {
            col_data$HIV <- ifelse(col_data$condition %in% c('HIV', 'HIV+OUD'), 1, 0)
            col_data$OUD <- ifelse(col_data$condition %in% c('OUD', 'HIV+OUD'), 1, 0)
            col_data$condition <- ifelse(col_data$condition == "HIV+OUD", "HIV_OUD", col_data$condition)
        }
    
        # Convert to correct variable types
        col_data$condition <- factor(col_data$condition)
        col_data$age <- as.integer(col_data$age)
        col_data$sex <- factor(col_data$sex)
        col_data$age_scaled <- scale(col_data$age, center = TRUE)
    
        # Create DESeq2 dataset
        if (condition == 'HIVOUD_Interaction') {
            dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~sex + age_scaled + HIV + OUD + HIV:OUD)
        } else {
            dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~sex + age_scaled + condition)
        }
    
        # Pre-filter genes with low expression
        keep <- rowSums((counts(dds) >= 3)) > 0
        dds <- dds[keep, ]
    
        female_associated_genes <- c("XIST")
        dds <- dds[!rownames(dds) %in% female_associated_genes, ]
    
        # Scaling factor for filtering
        genes_before_filter <- nrow(dds)
        scaling_factor <- 27000 / genes_before_filter
        min_threshold <- 70 / scaling_factor
    
        keep <- rowSums((counts(dds) >= min_threshold)) > 1
        dds <- dds[keep, ]
    
        genes_after_filter <- nrow(dds)
        print(paste("Genes before filtering:", genes_before_filter))
        print(paste("Genes after filtering:", genes_after_filter))
        print(paste("Minimum expression threshold used:", min_threshold))
    
        # Filter out donors if too few cells
        min_cell_count <- 10
        donor_counts <- colSums(counts(dds))
        dds <- dds[, donor_counts > min_cell_count]
    
        # PCA Analysis
        rld <- rlog(dds, blind=TRUE)
        ggsave(file.path(output_dir_2, paste(condition, cell_type, "pca_sex.png", sep = "_")), plot = plotPCA(rld, intgroup = "sex"))
        ggsave(file.path(output_dir_2, paste(condition, cell_type, "pca_age.png", sep = "_")), plot = plotPCA(rld, intgroup = "age"))
        ggsave(file.path(output_dir_2, paste(condition, cell_type,"pca_condition.png", sep = "_")), 
               plot = plotPCA(rld, intgroup = "condition"))
    
        # Heatmap
        rld_mat <- assay(rld)
        rld_cor <- cor(rld_mat)
        pheatmap(rld_cor, annotation_col = col_data["condition"], cluster_cols = FALSE, cluster_rows=FALSE)
        
        # Run DESeq2
        dds <- DESeq(dds)
    
        if (condition == 'HIV') {
            res <- lfcShrink(dds, coef = 'condition_HIV_vs_Control', type ='apeglm')
        } else if (condition == 'OUD') {
            res <- lfcShrink(dds, coef = 'condition_OUD_vs_Control', type ='apeglm')
        } else if (condition == 'HIVOUD_Interaction') {
            res <- lfcShrink(dds, coef = 'HIV1.OUD1', type ='apeglm')
        }
    
        res_tbl <- res %>%
        data.frame() %>%
        rownames_to_column(var = "gene") %>%
        as_tibble() %>%
        arrange(padj)
        
        write.csv(as.data.frame(res_tbl),
          file=paste(output_dir_2,"/",cell_type,"_",region,"_",condition,"_results_all_genes.csv", sep = ""))
        print("Saved results tables...")}
        
        # MA Plot
        png(file.path(output_dir_2, "/", cell_type, "_MA_plot.png"))
        plotMA(res, main = "DESeq2 MA Plot", ylim = c(-2, 2))
        dev.off()
        
        # Heatmap of top differentially expressed genes
        topGenes <- rownames(res)[1:20]
        norm_counts <- assay(vst(dds))[topGenes, ]
        
        png(file.path(output_dir_2,"/", cell_type, "_Heatmap.png"))
        pheatmap(norm_counts, cluster_rows = TRUE, cluster_cols = TRUE, 
                 color = colorRampPalette(brewer.pal(9, "Blues"))(50))
        dev.off()
    }

