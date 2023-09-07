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
########### Script Arguments
################################################################################
# Load annotation file
anno_file <- read_csv("Data/GRCm39_rel104_full_annotation.txt")

# Load count file
count_file <- "Data/raw_counts.csv"
count_matrix <- read_csv(count_file)

# Load column metadata
metadata_file <- "Data/column_data.xlsx"
coldata <- read_xlsx(metadata_file)
# Set column names to include as factor
colnames(coldata)
factor_levels <- c("condition","group","genotype","tissue","type")

# Set design formula
design_factors <- c("group")
design_formula <- as.formula(paste("~",design_factors))
print(design_formula)

# Set samples to remove
samples_to_remove <- c(23, 28) 

# Set pre-filter cutoff (counts)
prefilter_cutoff <- 10
# Set padj and log2foldchange cutoff to be a DE_gene
sign_cutoff <- 0.05
fc_cutoff <- 0

# Create comparison information
comp_grouping <- "group"


################################################################################
########### Load Count_Matrix and order it correctly
################################################################################
# Order columns according to there sample number
# store count data columns
data_matrix <- count_matrix[-1]
# store name and data the first column (Ensemble_id)
gene_col_name <- colnames(count_matrix)[1]
gene_col <- count_matrix[1]

# Get column names
columns <- colnames(count_matrix)[-1]
# Get there sample number
column_sn <- c()
for (i in str_split(columns,"_")){column_sn <- append(column_sn,as.numeric(i[1]))}
if (length(is.na(column_sn)) > 0){
  print(paste0("The sample number was not found! Maybe split it from the name using other split string"))
  column_sn <- c()
  for (i in str_split(columns,"-")){column_sn <- append(column_sn,as.numeric(i[1]))}
}
# Order these sample numbers
column_sn_ordered <- sort(column_sn)
# Fetch index position of these sorted numbers in pre-sorted sample number list
index <- c()
for (n in column_sn_ordered){index <- append(index,which(column_sn == n))}
# Apply to the columnames and add first column
ordered_count_mtx <- data_matrix[,index]
ordered_count_mtx <- add_column(gene_col,ordered_count_mtx)
# Convert to a data.frame
ordered_count_mtx <- as.data.frame(ordered_count_mtx)
ordered_count_mtx <- column_to_rownames(ordered_count_mtx,var = gene_col_name)
coldata$sample <- colnames(ordered_count_mtx)

# Inspect
view(ordered_count_mtx)

# Clean up
rm(data_matrix,count_matrix,gene_col,column_sn,column_sn_ordered,columns,i,n,index,gene_col_name)
gc()

################################################################################
########### Set factor levels of count_matrix metadata
################################################################################
# Setting factors
for (factor in factor_levels){
  coldata[[factor]] <- as.factor(coldata[[factor]])
}
# Clean up
rm(factor)
gc()



################################################################################
########### CREATE DESEQ2 OBJECT
################################################################################
# Check if coldata is the same as the count_matrix
if (all(coldata$sample %in% colnames(ordered_count_mtx))){
  # Check if samples are ordered the same
  if (all(coldata$sample == colnames(ordered_count_mtx))){
    # Create object
    print("coldata named and ordered in the same way as the ordered_count_mtx")
    dds <- DESeqDataSetFromMatrix(ordered_count_mtx,
                                  colData = coldata,
                                  design = design_formula) # MAKE CHANGEABLE IN PYTHON ~group
  } else {
    # Order samples depending on coldata
    print("The sample names are not in the same order and will be rearranged!")
    ordered_count_mtx <- ordered_count_mtx[, coldata$sample]
    # Create object
    dds <- DESeqDataSetFromMatrix(ordered_count_mtx,
                                  colData = coldata,
                                  design = design_formula) # MAKE CHANGEABLE IN PYTHON ~group
  }
  
} else {
  print("The coldata and ordered_count_mtx are different!")
  # show differences
  coldata$sample %in% colnames(ordered_count_mtx)
}



################################################################################
########### REMOVE UNWANTED SAMPLES
################################################################################
# Remove unwanted samples
if (length(samples_to_remove != 0 )){
  samples_to_remove <- as.integer(samples_to_remove)
  dds <- dds[, ! dds$sample_number %in% samples_to_remove]
  print(samples_to_remove)
  print(colData(dds))

  # remove old metadata for all factors
  for (columns in names(Filter(is.factor, colData(dds)))){
    dds[[columns]] <- droplevels(dds[[columns]])
  }
}

################################################################################
########### Quality Control
################################################################################
### Inspect samples ###
# Variance Stabilizing transformation
vsd <- vst(dds, blind = F)
# extract the vst matris from the object
vsd_mat <- assay(vsd)
# compute pairwise correlation values
vsd_cor <- cor(vsd_mat)

# Heatmap
pheatmap(vsd_cor)

# PCA
colnames(colData(dds))
plotPCA(vsd, intgroup = "group") + theme_classic()



################################################################################
########### PRE-FILTERING (cpm or count based)
################################################################################
# Get number of genes before filtering
num_genes_nofilter <- dim(dds)[1]
# Get lowest amount of occurrences of a condition
num_samples <- min(summary(dds[[comp_grouping]]))
# Keep genes that have atleast a cpm of x in the minimum amount of samples in a condition (better then just filter on the counts)
keep <- rowSums(counts(dds) >= prefilter_cutoff) >= num_samples
filtered_dds <- dds[keep,]
num_genes_filter <- dim(filtered_dds)[1]
print(paste0("Genes Before Filtering: ",num_genes_nofilter," | Genes After Filtering: ",num_genes_filter))

# Clean up
rm(num_genes_filter,num_genes_nofilter,keep)
gc()


################################################################################
########### Save DESeq2 Object
################################################################################
saveRDS(filtered_dds, "Data/DESeq2_object.rds")
