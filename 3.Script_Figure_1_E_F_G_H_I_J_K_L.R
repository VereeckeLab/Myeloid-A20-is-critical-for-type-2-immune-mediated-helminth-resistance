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
library(GSVA)

source("2.Script_Functions.R")

################################################################################
########### Load DESeq2 Object and needed arguments
################################################################################
# Load dds object created in 1.Script_DESeq2_Pre_processing.R
filtered_dds <- readRDS("Data/DESeq2_object.rds")
dds <- readRDS("Data/DESeq2_unfiltered_object.rds")

# Load annotation file
anno_file <- read_csv("Data/GRCm39_rel104_full_annotation.txt")

# File containing M1/M2 linked genes
data_in <- read_excel(path = "Data/M1_M2_genes.xlsx")
# MsigDB C2 genesets
MsigdbC2 <- readRDS("Data/Mm.c2.all.v7.1.entrez.rds")

# Set cutoffs
prefilter_cutoff <- 10
sign_cutoff <- 0.05
fc_cutoff <- 0

# define factors
factor_levels <- c("condition","group","genotype","tissue","type")


################################################################################
########### DE analysis colon A20 vs WT (group D vs B)
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

################################################################################
########### GSEA (Figure 1G)
################################################################################
#### GSEA
shrunk_results <- read_xlsx("Data/Results/colon_A20_vs_WT_GSEA_filtered.xlsx")
DE_shrunk_results <- shrunk_results[shrunk_results$DE_gene == T,]
# Generate ranked gene file
ranks_lfc <- shrunk_results[order(shrunk_results$log2FoldChange,decreasing = T),c("Entrez_GSEA","log2FoldChange")] # GSEA YOU CAN DO FROM NON DE GENES OR FROM DE GENES
ranks_lfc_DE <- DE_shrunk_results[order(DE_shrunk_results$log2FoldChange,decreasing = T),c("Entrez_GSEA","log2FoldChange")]
ranks_lfc_DE <- deframe(ranks_lfc_DE)
# If duplicate gene names present, average the values
if( sum(duplicated(ranks_lfc$Entrez_GSEA)) > 0) {
  ranks_lfc <- aggregate(.~Entrez_GSEA, FUN = mean, data = ranks_lfc)
  ranks_lfc <- ranks_lfc[order(ranks_lfc$log2FoldChange, decreasing = T),]
}
# Turn the dataframe into a named vector for fgsea()
ranks_lfc <- deframe(ranks_lfc)

# Control how specific the analysis is with minSize and maxSize
fgsea_results <- fgsea(MsigdbC2,
                       ranks_lfc,
                       minSize = 3,
                       maxSize = 1000)
# Full Figure
GSEA_plot(fgsea_results,npathw = 40)
# Select pathways from the top significant results
test_set <- fgsea_results[fgsea_results$pathway %in% c(#"REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
                                                       "FOSTER_TOLERANT_MACROPHAGE_UP",
                                                       "BROWNE_INTERFERON_RESPONSIVE_GENES",
                                                       "PID_IL12_2PATHWAY",
                                                       #"ZHANG_INTERFERON_RESPONSE",
                                                       "BOSCO_TH1_CYTOTOXIC_MODULE",
                                                       "SANA_TNF_SIGNALING_UP",
                                                       "PID_IL23_PATHWAY",
                                                       "SANA_RESPONSE_TO_IFNG_UP",
                                                       #"GRANDVAUX_IRF3_TARGETS_UP",
                                                       #"REACTOME_IRF3_MEDIATED_INDUCTION_OF_TYPE_I_IFN",
                                                       "SEMENZA_HIF1_TARGETS",
                                                       #"BIOCARTA_NO2IL12_PATHWAY",
                                                       "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
                                                       "REACTOME_INFLAMMASOMES",
                                                       "TAVOR_CEBPA_TARGETS_UP",
                                                       "REACTOME_INTERFERON_GAMMA_SIGNALING",
                                                       #"REACTOME_INTERLEUKIN_1_FAMILY_SIGNALING",
                                                       #"REACTOME_TOLL_LIKE_RECEPTOR_CASCADES",
                                                       "PID_IL27_PATHWAY",
                                                       #"MARTIN_NFKB_TARGETS_UP",
                                                       "COATES_MACROPHAGE_M1_VS_M2_DN",
                                                       #"REACTOME_TOLL_LIKE_RECEPTOR_TLR1_TLR2_CASCADE",
                                                       #"REACTOME_TOLL_LIKE_RECEPTOR_4_TLR4_CASCADE",
                                                       "REACTOME_INTERLEUKIN_1_SIGNALING",
                                                       "MARTINELLI_IMMATURE_NEUTROPHIL_DN"
                                                       #"REACTOME_TRAF6_MEDIATED_NF_KB_ACTIVATION",
                                                       #"AUJLA_IL22_AND_IL17A_SIGNALING"
                                                       )]
p3 <- GSEA_plot(test_set,npathw = 10000,NES_cutoff = 0,color_down = "#ff1493",title = "")
ggsave("Data/Figures/colon_A20_vs_WT_GSEA.png",plot = p3,units = "px", width = 900, height = 1200, dpi = 85)

################################################################################
########### GSVA DATA
################################################################################
# fetch normalized counts
dds <- estimateSizeFactors(dds)
norm.counts <- as.data.frame(counts(dds,normalized = T))
# add annotation
annotation <- anno_file
norm.counts <- rownames_to_column(norm.counts,var = "Gene.stable.ID")
norm.counts <- merge(x = norm.counts, y = annotation[,c("Gene.stable.ID","Gene.name")], by = "Gene.stable.ID" , all.x = T)
# remove gene counts
norm.counts$Gene.stable.ID <- NULL
# remove empty gene.name columns
table(norm.counts$Gene.name != "")
norm.counts <- norm.counts[norm.counts$Gene.name != "",]
# remove na
table(is.na(norm.counts$Gene.name))
norm.counts <- norm.counts[!is.na(norm.counts$Gene.name),]

# average out duplicates
if( sum(duplicated(norm.counts$Gene.name)) > 0) {
  norm.counts <- aggregate(.~Gene.name, FUN = mean, data = norm.counts)
}
# add gene names to rownames
norm.counts <- column_to_rownames(norm.counts,var = "Gene.name")
# change to matrix
norm.counts <- as.matrix(norm.counts)

# fetch gene sets M1 and M2
genes_m1 <- as.vector(na.omit(data_in$M1))
genes_m2 <- as.vector(na.omit(data_in$M2))
geneset <- list(genes_m1,genes_m2)
names(geneset) <- c("CLASSICAL ACTIVATION OF MACROPHAGES","ALTERNATIVE ACTIVATION OF MACROPHAGES")

# GSVA
gsva_results <- gsva(norm.counts,geneset, method = "gsva")
score <- as.data.frame(t(gsva_results))
score <- rownames_to_column(score, var = "samples")
# add metadata
score$Condition <- colData(dds)$condition
levels(score$Condition)
score$Condition <-  factor(score$Condition,levels = c("colon_steady_state_wt","colon_steady_state_A20","colon_trichuris_wt","colon_trichuris_A20","joint_steady_state_wt","joint_steady_state_A20"))
score$Genotype <- colData(dds)$genotype 
levels(score$Genotype) <- c("A20\u1d50\u02b8\u1d49\u02e1\u207b\u1d37\u1d3c","WT") #"\u0394-A20"
score$Genotype <-  factor(score$Genotype,levels = c("WT","A20\u1d50\u02b8\u1d49\u02e1\u207b\u1d37\u1d3c"))
score$name <- score$Condition
levels(score$name) <- c("Colon Steady State","Colon Steady State","Colon T. muris","Colon T. muris","Joint Steady State","Joint Steady State")
score$Tissue <- score$Condition
levels(score$Tissue) <- c("Colon","Colon","Colon","Colon","Joint","Joint")
score$State <- score$Condition
levels(score$State) <- c("Steady State","Steady State","T. muris","T. muris","Steady State","Steady State")
score$Shape <- score$Condition
levels(score$State) <- c("dot","circle","dot","dot","dot","circle")
score$Color <- score$Condition
levels(score$State) <- c("black","black","blue","red","black","black")

################################################################################
########### GSVA (Figure 1H)
################################################################################
# ONLY SS COLON M1
score.colon <- score[score$Condition %in% c("colon_steady_state_wt","colon_steady_state_A20"),]
ggplot(score.colon,aes(x= Genotype, y= .data[["CLASSICAL ACTIVATION OF MACROPHAGES"]], group = Genotype , fill = name, shape = Shape )) +
  stat_boxplot(geom ='errorbar',width = 0.25) + # horizontal lines on boxplot print first
  geom_boxplot(outlier.shape = NA, width = 0.6, lwd = 0.9) + facet_grid(. ~ name) +
  scale_x_discrete(expand=c(0.8,0)) +
  scale_shape_manual(values=c(16,1,16,16,16,1) , guide = "none") +
  scale_fill_manual(values=c("white","white", "white")) +
  ggtitle("CLASSICAL ACTIVATION") +
  ylab("GSVA Enrichment Score") + xlab("") +
  ylim(c(-0.45,0.21)) +
  geom_point(aes(colour = Color), position=position_jitterdodge(dodge.width = 0.75), size = 6.5) +
  scale_colour_manual(values=c("black", "black","blue", "red","black", "black"), guide = "none")+
  theme_classic() +
  theme(text = element_text(size = 26, color = "black"),
        axis.text = element_text(size = 21, color = "black"),
        axis.text.x = element_text(size = 32, color = "black"),
        axis.title.y = element_text(size = 25, color = "black"),
        strip.text.x = element_text(size = 25),
        strip.background = element_rect(fill="grey"),
        #legend.title = element_text(face = "bold"),
        plot.title = element_text(size = 27,hjust = 0.5)) +
  guides(fill = "none") +
  stat_compare_means(comparisons = list(c("WT", "A20\u1d50\u02b8\u1d49\u02e1\u207b\u1d37\u1d3c")), label = "p.signif", size = 18, label.y = 0.14,vjust = -0.3, bracket.size = 0.9) +
  stat_compare_means(comparisons = list(c("WT", "A20\u1d50\u02b8\u1d49\u02e1\u207b\u1d37\u1d3c")), size = 8, label.y = 0.14)
ggsave("Data/Figures/colon_A20_vs_WT_M1_SS_GSVA.png", units = "px", width = 500, height = 1100, dpi=85)


################################################################################
########### DE analysis Joint A20 vs WT (group D vs B)
################################################################################
# Subset dds object
subset_dds <- filtered_dds[,filtered_dds[["group"]] %in% c("C","A")] 
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
subset_dds[["group"]] <- relevel(subset_dds[["group"]],ref = "A")
# Perform pre-processing
f_dds_DE <- DESeq(subset_dds)
# Get expression results
DE_result <- results(f_dds_DE,
                     contrast=c("group","C","A"),
                     alpha = sign_cutoff,
                     lfcThreshold = fc_cutoff)
summary(DE_result)
# Save results
save_all_DE_results(dds_DE = f_dds_DE,
                    dds_res = DE_result,
                    save_dir = "Data/Results",
                    annotation = anno_file,
                    name = "joint_A20_vs_WT")

# Fold change shrinkage
# Shrink DE results -> adjusted LogFC depending on sample variation
# DESeq2 shrinks the log2 fold change (LFC) values of gene expression changes across conditions
# to address the problem of overestimation of fold changes in RNA-seq experiments with low sample sizes or low read counts.
DE_result_shrunk <- lfcShrink(dds = f_dds_DE,
                                contrast=c("group","C","A"),
                                type="ashr",
                                res = DE_result)
summary(DE_result_shrunk)

save_all_DE_results(dds_DE = f_dds_DE,
                    dds_res = DE_result_shrunk,
                    save_dir = "Data/Results",
                    annotation = anno_file,
                    name_file = "joint_A20_vs_WT_GSEA")


################################################################################
########### VolcanoPlot (Figure 1I)
################################################################################
# Load DE results
results <- read_xlsx("Data/Results/joint_A20_vs_WT_filtered.xlsx")
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
                       x_scale_breaks = c(-8,-5,-2,0,2,5,8),
                       x_axis_limits = c(-10,10),
                       logFC_cutoff = 0,
                       point_size = 3,
                       text_title_size = 30)
ggsave("Data/Figures/joint_A20_vs_WT.png",units = "px", width = 960, height = 1200, dpi = 85)


################################################################################
########### Heatmap (Figure 1J)
################################################################################
# M1 and M2
dds_vst <- vst(f_dds_DE, blind = TRUE)
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
keep <- na.omit(coldata[coldata$tissue == "JOINT", "sample" ])
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
annoCol<-list(condition=c(joint_steady_state_A20="red", joint_steady_state_wt="blue"))

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
ggsave("Data/Figures/joint_A20_vs_WT_M1_heatmap.png",plot = p,units = "px", width = 650, height = 1200, dpi = 85)


################################################################################
########### GSEA (Figure 1K)
################################################################################
#### GSEA
shrunk_results <- read_xlsx("Data/Results/joint_A20_vs_WT_GSEA_filtered.xlsx")
DE_shrunk_results <- shrunk_results[shrunk_results$DE_gene == T,]
# Generate ranked gene file
ranks_lfc <- shrunk_results[order(shrunk_results$log2FoldChange,decreasing = T),c("Entrez_GSEA","log2FoldChange")] # GSEA YOU CAN DO FROM NON DE GENES OR FROM DE GENES
ranks_lfc_DE <- DE_shrunk_results[order(DE_shrunk_results$log2FoldChange,decreasing = T),c("Entrez_GSEA","log2FoldChange")]
ranks_lfc_DE <- deframe(ranks_lfc_DE)
# If duplicate gene names present, average the values
if( sum(duplicated(ranks_lfc$Entrez_GSEA)) > 0) {
  ranks_lfc <- aggregate(.~Entrez_GSEA, FUN = mean, data = ranks_lfc)
  ranks_lfc <- ranks_lfc[order(ranks_lfc$log2FoldChange, decreasing = T),]
}
# Turn the dataframe into a named vector for fgsea()
ranks_lfc <- deframe(ranks_lfc)

# Control how specific the analysis is with minSize and maxSize
fgsea_results <- fgsea(MsigdbC2,
                       ranks_lfc,
                       minSize = 3,
                       maxSize = 1000)
# Full Figure
GSEA_plot(fgsea_results,npathw = 40)
# Select pathways from the top significant results
test_set <- fgsea_results[fgsea_results$pathway %in% c("DASU_IL6_SIGNALING_UP",
                                                       #"REACTOME_INTERLEUKIN_12_SIGNALING",
                                                       "AUJLA_IL22_AND_IL17A_SIGNALING",
                                                       #"BIOCARTA_STAT3_PATHWAY",
                                                       #"BASSO_CD40_SIGNALING_UP",
                                                       #"REACTOME_INTERLEUKIN_21_SIGNALING",
                                                       #"WANG_TNF_TARGETS",
                                                       #"BIOCARTA_NOS1_PATHWAY",
                                                       "REACTOME_INFLAMMASOMES",
                                                       #"PID_IL27_PATHWAY",
                                                       #"REACTOME_TRAF6_MEDIATED_IRF7_ACTIVATION",
                                                       #"REACTOME_TOLL_LIKE_RECEPTOR_CASCADES",
                                                       #"REACTOME_INTERLEUKIN_1_FAMILY_SIGNALING",
                                                       "REACTOME_INTERLEUKIN_1_SIGNALING",
                                                       "BROWN_MYELOID_CELL_DEVELOPMENT_UP",
                                                       #"FOSTER_TOLERANT_MACROPHAGE_DN",
                                                       #"CROONQUIST_IL6_DEPRIVATION_DN",
                                                       #"PID_TOLL_ENDOGENOUS_PATHWAY",
                                                       "NEMETH_INFLAMMATORY_RESPONSE_LPS_UP",
                                                       "BIOCARTA_NO2IL12_PATHWAY",
                                                       "MARTIN_NFKB_TARGETS_UP",
                                                       "ZHENG_IL22_SIGNALING_UP",
                                                       #"REACTOME_TRAF6_MEDIATED_NF_KB_ACTIVATION",
                                                       #"BOSCO_TH1_CYTOTOXIC_MODULE",
                                                       #"DEBOSSCHER_NFKB_TARGETS_REPRESSED_BY_GLUCOCORTICOIDS",
                                                       "PID_IL23_PATHWAY",
                                                       "TAVOR_CEBPA_TARGETS_UP",
                                                       "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
                                                       "PID_IL12_2PATHWAY",
                                                       #"KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
                                                       #"REACTOME_INTERLEUKIN_10_SIGNALING",
                                                       #"SANA_TNF_SIGNALING_UP",
                                                       "ZHANG_RESPONSE_TO_IKK_INHIBITOR_AND_TNF_UP",
                                                       #"SANA_RESPONSE_TO_IFNG_UP",
                                                       #"LINDSTEDT_DENDRITIC_CELL_MATURATION_A",
                                                       #"MOSERLE_IFNA_RESPONSE",
                                                       "REACTOME_INTERFERON_SIGNALING",
                                                       "REACTOME_INTERFERON_GAMMA_SIGNALING"
                                                       #"SEKI_INFLAMMATORY_RESPONSE_LPS_UP"
                                                       )]
p3 <- GSEA_plot(test_set,npathw = 10000,NES_cutoff = 0,color_down = "#ff1493",title = "")
ggsave("Data/Figures/joint_A20_vs_WT_GSEA.png",plot = p3,units = "px", width = 900, height = 1200, dpi = 85)


################################################################################
########### GSVA (Figure 1L)
################################################################################
# ONLY SS JOINT M1
score.joint <- score[score$Condition %in% c("joint_steady_state_wt","joint_steady_state_A20"),]
ggplot(score.joint,aes(x= Genotype, y= .data[["CLASSICAL ACTIVATION OF MACROPHAGES"]], group = Genotype , fill = name, shape = Shape )) +
  stat_boxplot(geom ='errorbar',width = 0.25) + # horizontal lines on boxplot print first
  geom_boxplot(outlier.shape = NA, width = 0.6, lwd = 0.9) + facet_grid(. ~ name) +
  scale_x_discrete(expand=c(0.8,0)) +
  scale_shape_manual(values=c(16,1,16,16,16,1) , guide = "none") +
  scale_fill_manual(values=c("white","white", "white")) +
  ggtitle("CLASSICAL ACTIVATION") +
  ylab("GSVA Enrichment Score") + xlab("") +
  ylim(c(-0.35,0.49)) + # Influences stats if a dot is not plotted!
  geom_point(aes(colour = Color), position=position_jitterdodge(dodge.width = 0.75), size = 6.5) +
  scale_colour_manual(values=c("black", "black","blue", "red","black", "black"), guide = "none")+
  theme_classic() +
  theme(text = element_text(size = 26, color = "black"),
        axis.text = element_text(size = 21, color = "black"),
        axis.text.x = element_text(size = 32, color = "black"),
        axis.title.y = element_text(size = 25, color = "black"),
        strip.text.x = element_text(size = 25),
        strip.background = element_rect(fill="grey"),
        #legend.title = element_text(face = "bold"),
        plot.title = element_text(size = 27, hjust = 0.5)) +
  guides(fill = "none") +
  stat_compare_means(comparisons = list(c("WT", "A20\u1d50\u02b8\u1d49\u02e1\u207b\u1d37\u1d3c")), label = "p.signif", size = 18, label.y = 0.4,vjust = -0.3, bracket.size = 0.9) +
  stat_compare_means(comparisons = list(c("WT", "A20\u1d50\u02b8\u1d49\u02e1\u207b\u1d37\u1d3c")), size = 8, label.y = 0.4)
ggsave("Data/Figures/joint_A20_vs_WT_M1_SS_GSVA.png", units = "px", width = 500, height = 1100, dpi=85)

# Validate statistics
wilcox.test(score.joint[score.joint$Genotype == "WT","CLASSICAL ACTIVATION OF MACROPHAGES"],
            score.joint[score.joint$Genotype == "A20\u1d50\u02b8\u1d49\u02e1\u207b\u1d37\u1d3c","CLASSICAL ACTIVATION OF MACROPHAGES"])

wilcox.test(score.joint[score.joint$Genotype == "WT","ALTERNATIVE ACTIVATION OF MACROPHAGES"],
            score.joint[score.joint$Genotype == "A20\u1d50\u02b8\u1d49\u02e1\u207b\u1d37\u1d3c","ALTERNATIVE ACTIVATION OF MACROPHAGES"])

