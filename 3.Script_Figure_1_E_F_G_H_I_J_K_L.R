################################################################################
########### Load Packages
################################################################################
library(readxl)
library(writexl)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(gtools)
library(pheatmap)
library(RColorBrewer)
library(grDevices)
library(ggpubr)
library(ggforce)
library(ggrepel)
library(DESeq2)
library(edgeR)
library(apeglm)
library(ashr)
library(org.Mm.eg.db)
library(fgsea)
library(gplots)

source("2.Script_Functions.R")

################################################################################
########### Load DESeq2 Object and needed arguments
################################################################################
# Load dds object created in 1.Script_DESeq2_Pre_processing.R
filtered_dds <- readRDS("Data/DESeq2_object.rds")

# Load annotation file
anno_file <- read_csv("Data/GRCm39_rel104_full_annotation.txt")

# Set cutoffs
prefilter_cutoff <- 10
sign_cutoff <- 0.05
fc_cutoff <- 0

# define factors
factor_levels <- c("condition","group","genotype","tissue","type")


################################################################################
########### DE analysis colon A20 vs WT (group A vs E)
################################################################################
# Subset dds object
subset_dds <- filtered_dds[,filtered_dds[["group"]] %in% c("D","B")] 
# Drop unused levels to avoid errors
for (factor in factor_levels){
  subset_dds[[factor]] <- droplevels(subset_dds[[factor]])
}
# quick check
head(subset_dds@colData)
vsd <- vst(subset_dds, blind = F)
# QC
plotPCA(vsd, intgroup = "condition") + theme_classic()

# Find differentially expressed genes
# Set reference group
subset_dds[["group"]] <- relevel(subset_dds[["group"]],ref = "B")
# Perform pre-processing
f_dds_DE <- DESeq(subset_dds)
# Get expression results
DE_result <- results(f_dds_DE,
                     contrast=c("group","D","B"),
                     alpha = sign_cutoff,
                     lfcThreshold = fc_cutoff)
summary(DE_result)
# Save results
save_all_DE_results(dds_DE = f_dds_DE,
                    dds_res = DE_result,
                    save_dir = "Data/Results",
                    annotation = anno_file,
                    name = "colon_A20_vs_WT")

# Fold change shrinkage
# Shrink DE results -> adjusted LogFC depending on sample variation
# DESeq2 shrinks the log2 fold change (LFC) values of gene expression changes across conditions
# to address the problem of overestimation of fold changes in RNA-seq experiments with low sample sizes or low read counts.
DE_result_shrunk <- lfcShrink(dds = f_dds_DE,
                                contrast=c("group","D","B"),
                                type="ashr",
                                res = DE_result)
summary(DE_result_shrunk)

save_all_DE_results(dds_DE = f_dds_DE,
                    dds_res = DE_result_shrunk,
                    save_dir = "Data/Results",
                    annotation = anno_file,
                    name_file = "colon_A20_vs_WT_GSEA")

################################################################################
########### VolcanoPlot (Figure 1E)
################################################################################
# Load DE results
results <- read_xlsx("Data/Results/colon_A20_vs_WT_filtered.xlsx")
# Remove low base mean? (extra filtering)
results <- results[results$baseMean > 10,]
# Number of up-regulated genes
sign_up <- results[results$DE_gene==T & results$log2FoldChange > 0,]
# Number of down-regulated genes
sign_down <- results[results$DE_gene==T & results$log2FoldChange < 0,]

top10_up <- head(sign_up[order(sign_up$log2FoldChange,decreasing = T),"Gene.name"],10)
top10_down <- head(sign_down[order(sign_down$log2FoldChange,decreasing = F),"Gene.name"],10)
top20 <- c(deframe(top10_up),deframe(top10_down))

get_paper_volcano_plot(results,
                       title = "",
                       important_genes = top20,
                       crop_value_cutoff = 1e-15,
                       colors = c("#00adef", "black", "#ff1493"),
                       x_scale_breaks = c(-6,-4,-2,0,2,4,6),
                       x_axis_limits = c(-10,10),
                       logFC_cutoff = 0,
                       point_size = 3,
                       text_title_size = 30)
ggsave("Data/Figures/colon_A20_vs_WT.png",units = "px", width = 960, height = 1200, dpi = 85)

################################################################################
########### Heatmap (Figure 1F)
################################################################################
# M1 and M2
dds_vst <- vst(f_dds_DE, blind = TRUE)
data_in <- read_excel(path = "Data/M1_M2_genes.xlsx")
genes_m1 <- as.vector(na.omit(data_in$M1))
genes_m2 <- as.vector(na.omit(data_in$M2))
# Expand rowData(dds_vst) with ensembl_id and gene_name
rowData(dds_vst)[["ensembl_id"]] <- rownames(rowData(dds_vst))
rowData(dds_vst)[["gene_name"]] <- anno_file$Gene.name[match(rownames(rowData(dds_vst)), anno_file$Gene.stable.ID)]

DE_genes_M1 <- deframe(results[(results$Gene.name %in% genes_m1 & results$DE_gene == T), "Gene.name" ])
DE_genes_M2 <- deframe(results[(results$Gene.name %in% genes_m2 & results$DE_gene == T), "Gene.name" ])

# Classically activated macrophages (M1) 
M1_dds_vst_DE <- dds_vst[rowData(dds_vst)[["gene_name"]] %in% DE_genes_M1, ] # important genes genes_m1 contains gene names
rownames(M1_dds_vst_DE) <- rowData(M1_dds_vst_DE)[["gene_name"]] #set rownames either "ensembl_id" or "gene_name"
# Alternatively activated macrophages (M2)
M2_dds_vst_DE <- dds_vst[rowData(dds_vst)[["gene_name"]] %in% DE_genes_M2, ] # important genes genes_m2 contains gene names
rownames(M2_dds_vst_DE) <- rowData(M2_dds_vst_DE)[["gene_name"]] #set rownames either "ensembl_id" or "gene_name"

# Get top M1 genes based on DEG results
mat <- assay(M1_dds_vst_DE)
head(mat)
head(mat)
coldata <- colData(M1_dds_vst_DE)
colnames(mat) <- gsub(pattern = "^X", replacement = "", x = colnames(mat))
colnames(mat) <- coldata[,"sample"]
keep <- na.omit(coldata[coldata$tissue == "COLON", "sample" ])
mat <- rownames_to_column(as.data.frame(mat), "Gene.name")
top <- mat
colnames(top)
top <- top[,c(keep,"Gene.name")]
top <- column_to_rownames(top, "Gene.name")
top <- top[!(rownames(top) %in% c("Cxcl3","H2-Ob","Acpp","Smim3","Dusp1","Tmem119","H2.Ob")),]

# create heatmap anno df
anno <- as.data.frame(colData(M1_dds_vst_DE)[, c("condition","sample"), drop = FALSE])
rownames(anno) <- NULL
anno <- column_to_rownames(anno, "sample")
levels(colData(M1_dds_vst_DE)$condition)
annoCol<-list(condition=c(colon_steady_state_A20="red", colon_steady_state_wt="blue"))

# set collor palettes
pal <- colorpanel(256, "#00adef", "white", "#ff1493")


# order heatmap according to there log2FC in the results (top)
datas <- assay(M1_dds_vst_DE)
datas <- datas[!(rownames(datas) %in% c("Cxcl3","H2-Ob","Acpp","Smim3","Dusp1","Tmem119","H2.Ob")),]
res_genes <- results[results$Gene.name %in% rownames(datas),]
nrow(res_genes) == nrow(top) # needs to be TRUE
res_genes <- res_genes[order(res_genes$log2FoldChange, decreasing = TRUE),"Gene.name"]
#res_genes <- make.names(res_genes) 
top_most_var <- top[match(deframe(res_genes), rownames(top)),]

p <- pheatmap(top_most_var[1:25,],
         cluster_cols=F,
         color = pal,
         #annotation_col= anno, #add grouping above heatmap
         #annotation_colors = annoCol,
         cluster_rows = F,
         show_rownames = T,
         show_colnames = F,
         #border_color = NA,
         fontsize = 16,
         scale = "row",
         fontsize_row = 20,
         height = 20,
         main = "")
ggsave("Data/Figures/colon_A20_vs_WT_M1_heatmap.png",plot = p,units = "px", width = 650, height = 1200, dpi = 85)
