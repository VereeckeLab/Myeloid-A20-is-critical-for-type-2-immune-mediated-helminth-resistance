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
########### DE analysis colon Trichuris muris A20 vs WT (group F vs E)
################################################################################
# Subset dds object
subset_dds <- filtered_dds[,filtered_dds[["group"]] %in% c("F","E")] 
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
subset_dds[["group"]] <- relevel(subset_dds[["group"]],ref = "E")
# Perform pre-processing
f_dds_DE <- DESeq(subset_dds)
# Get expression results
DE_result <- results(f_dds_DE,
                     contrast=c("group","F","E"),
                     alpha = sign_cutoff,
                     lfcThreshold = fc_cutoff)
summary(DE_result)
# Save results
save_all_DE_results(dds_DE = f_dds_DE,
                    dds_res = DE_result,
                    save_dir = "Data/Results",
                    annotation = anno_file,
                    name = "colon_TM_A20_vs_WT")

# Fold change shrinkage
# Shrink DE results -> adjusted LogFC depending on sample variation
# DESeq2 shrinks the log2 fold change (LFC) values of gene expression changes across conditions
# to address the problem of overestimation of fold changes in RNA-seq experiments with low sample sizes or low read counts.
DE_result_shrunk <- lfcShrink(dds = f_dds_DE,
                                contrast=c("group","F","E"),
                                type="ashr",
                                res = DE_result)
summary(DE_result_shrunk)

save_all_DE_results(dds_DE = f_dds_DE,
                    dds_res = DE_result_shrunk,
                    save_dir = "Data/Results",
                    annotation = anno_file,
                    name_file = "colon_TM_A20_vs_WT_GSEA")


################################################################################
########### VolcanoPlot (Figure 4A)
################################################################################
# Load DE results
results <- read_xlsx("Data/Results/colon_TM_A20_vs_WT_filtered.xlsx")
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
                       x_scale_breaks = c(-8,-6,-4,-2,0,2,4,6,8),
                       x_axis_limits = c(-10,10),
                       logFC_cutoff = 0, #0.58
                       point_size = 3,
                       text_title_size = 30)
ggsave("Data/Figures/colon_TM_A20_vs_WT.png",units = "px", width = 960, height = 1200, dpi = 85)


################################################################################
########### GSEA (Figure 4B)
################################################################################
#### GSEA
shrunk_results <- read_xlsx("Data/Results/colon_TM_A20_vs_WT_GSEA_filtered.xlsx")
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
test_set <- fgsea_results[fgsea_results$pathway %in% c(#"REACTOME_REGULATION_OF_BETA_CELL_DEVELOPMENT",
                                                       "PID_IL4_2PATHWAY",
                                                       #"CROONQUIST_IL6_DEPRIVATION_DN",
                                                       "LINDSTEDT_DENDRITIC_CELL_MATURATION_C",
                                                       "KEGG_INTESTINAL_IMMUNE_NETWORK_FOR_IGA_PRODUCTION",
                                                       "BIOCARTA_GATA3_PATHWAY",
                                                       "ZHANG_RESPONSE_TO_IKK_INHIBITOR_AND_TNF_DN", #MAKE THIS POSITIVE
                                                       #"HOFFMANN_IMMATURE_TO_MATURE_B_LYMPHOCYTE_DN",
                                                       "NAKAJIMA_EOSINOPHIL",
                                                       "KARAKAS_TGFB1_SIGNALING",
                                                       "REACTOME_ACTIVATION_OF_PPARGC1A_PGC_1ALPHA_BY_PHOSPHORYLATION",
                                                       "SHAFFER_IRF4_TARGETS_IN_ACTIVATED_B_LYMPHOCYTE",
                                                       #"WAKABAYASHI_ADIPOGENESIS_PPARG_BOUND_8D",
                                                       #"BAKER_HEMATOPOESIS_STAT1_TARGETS",
                                                       #"KEGG_JAK_STAT_SIGNALING_PATHWAY",
                                                       #"WIERENGA_STAT5A_TARGETS_DN",
                                                       "TIMOFEEVA_GROWTH_STRESS_VIA_STAT1_DN",
                                                       #"REACTOME_TRAF6_MEDIATED_INDUCTION_OF_TAK1_COMPLEX_WITHIN_TLR4_COMPLEX",
                                                       "WIERENGA_STAT5A_TARGETS_GROUP2",
                                                       "AUJLA_IL22_AND_IL17A_SIGNALING",
                                                       #"REACTOME_TRAF6_MEDIATED_IRF7_ACTIVATION",
                                                       #"REACTOME_MYD88_INDEPENDENT_TLR4_CASCADE",
                                                       "REACTOME_TOLL_LIKE_RECEPTOR_4_TLR4_CASCADE",
                                                       "LIAN_NEUTROPHIL_GRANULE_CONSTITUENTS",
                                                       "BIOCARTA_NO2IL12_PATHWAY",
                                                       "MARTINELLI_IMMATURE_NEUTROPHIL_DN",
                                                       "REACTOME_INTERLEUKIN_1_SIGNALING",
                                                       "REACTOME_IRF3_MEDIATED_INDUCTION_OF_TYPE_I_IFN",
                                                       "DER_IFN_GAMMA_RESPONSE_UP",
                                                       "COATES_MACROPHAGE_M1_VS_M2_DN",
                                                       #"REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
                                                       #"GALINDO_IMMUNE_RESPONSE_TO_ENTEROTOXIN",
                                                       "REACTOME_INFLAMMASOMES",
                                                       "SANA_TNF_SIGNALING_UP",
                                                       #"NEMETH_INFLAMMATORY_RESPONSE_LPS_UP",
                                                       "REACTOME_INNATE_IMMUNE_SYSTEM",
                                                       #"BROWNE_INTERFERON_RESPONSIVE_GENES",
                                                       "ZHANG_INTERFERON_RESPONSE",
                                                       "REACTOME_NEUTROPHIL_DEGRANULATION",
                                                       "FOSTER_TOLERANT_MACROPHAGE_UP",
                                                       "BROWN_MYELOID_CELL_DEVELOPMENT_UP"
                                                       #"ZHENG_IL22_SIGNALING_UP"
                                                       )]
p3 <- GSEA_plot(test_set,npathw = 10000,NES_cutoff = 0,title = "")
ggsave("Data/Figures/colon_TM_A20_vs_WT_GSEA.png",plot = p3,units = "px", width = 1200, height = 1200, dpi = 85)


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
# ONLY TM COLON M1
score.colon.tm <- score[score$Condition %in% c("colon_trichuris_wt","colon_trichuris_A20"),]
p1 <- ggplot(score.colon.tm,aes(x= Genotype, y= .data[["CLASSICAL ACTIVATION OF MACROPHAGES"]], group = Genotype , fill = name, shape = Shape )) +
  stat_boxplot(geom ='errorbar',width = 0.25) + # horizontal lines on boxplot print first
  geom_boxplot(outlier.shape = NA, width = 0.6, lwd = 0.9) + facet_grid(. ~ name) +
  scale_x_discrete(expand=c(0.8,0)) +
  scale_shape_manual(values=c(16,16) , guide = "none") +
  scale_fill_manual(values=c("white","white", "white")) +
  ggtitle("CLASSICAL ACTIVATION") +
  ylab("GSVA Enrichment Score") + xlab("") +
  ylim(c(-0.33,0.65)) + # Influences stats if a dot is not plotted!
  geom_point(aes(colour = Color), position=position_jitterdodge(dodge.width = 0.75), size = 6.5) +
  scale_colour_manual(values=c("blue", "red"), guide = "none")+
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
  stat_compare_means(comparisons = list(c("WT", "A20\u1d50\u02b8\u1d49\u02e1\u207b\u1d37\u1d3c")), label = "p.signif", size = 18, label.y = 0.48,vjust = -0.3, bracket.size = 0.9) +
  stat_compare_means(comparisons = list(c("WT", "A20\u1d50\u02b8\u1d49\u02e1\u207b\u1d37\u1d3c")), size = 8, label.y = 0.48)
p1
ggsave("Data/Figures/colon_A20_vs_WT_M1_TM_GSVA.png", units = "px", width = 750, height = 1100, dpi=85)

p2 <- ggplot(score.colon.tm,aes(x= Genotype, y= .data[["ALTERNATIVE ACTIVATION OF MACROPHAGES"]], group = Genotype , fill = name, shape = Shape )) +
  stat_boxplot(geom ='errorbar',width = 0.25) + # horizontal lines on boxplot print first
  geom_boxplot(outlier.shape = NA, width = 0.6, lwd = 0.9) + facet_grid(. ~ name) +
  scale_x_discrete(expand=c(0.8,0)) +
  scale_shape_manual(values=c(16,16) , guide = "none") +
  scale_fill_manual(values=c("white","white", "white")) +
  ggtitle("ALTERNATIVE ACTIVATION") +
  ylab("GSVA Enrichment Score") + xlab("") +
  ylim(c(-0.33,0.65)) + # Influences stats if a dot is not plotted!
  geom_point(aes(colour = Color), position=position_jitterdodge(dodge.width = 0.75), size = 6.5) +
  scale_colour_manual(values=c("blue", "red"), guide = "none")+
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
  stat_compare_means(comparisons = list(c("WT", "A20\u1d50\u02b8\u1d49\u02e1\u207b\u1d37\u1d3c")), label = "p.signif", size = 18, label.y = 0.48,vjust = -0.3, bracket.size = 0.9) +
  stat_compare_means(comparisons = list(c("WT", "A20\u1d50\u02b8\u1d49\u02e1\u207b\u1d37\u1d3c")), size = 8, label.y = 0.48)
p2
ggsave("Data/Figures/colon_A20_vs_WT_M2_TM_GSVA.png", units = "px", width = 750, height = 1100, dpi=85)

ggarrange(p1,p2,ncol = 1)
ggsave("Data/Figures/colon_A20_vs_WT_M1_M2_Combined_TM_GSVA.png", units = "px", width = 500, height = 1100, dpi=85)

# Check statistics
pt1 <- p1 + stat_compare_means(method = "wilcox.test")
pt2 <- p2 + stat_compare_means(method = "wilcox.test")
wilcox.test(score.colon.tm[score.colon.tm$Genotype == "WT","CLASSICAL ACTIVATION OF MACROPHAGES"],
            score.colon.tm[score.colon.tm$Genotype == "A20\u1d50\u02b8\u1d49\u02e1\u207b\u1d37\u1d3c","CLASSICAL ACTIVATION OF MACROPHAGES"])

wilcox.test(score.colon.tm[score.colon.tm$Genotype == "WT","ALTERNATIVE ACTIVATION OF MACROPHAGES"],
            score.colon.tm[score.colon.tm$Genotype == "A20\u1d50\u02b8\u1d49\u02e1\u207b\u1d37\u1d3c","ALTERNATIVE ACTIVATION OF MACROPHAGES"])
