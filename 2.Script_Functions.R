# Store DE results
save_all_DE_results <- function(dds_DE , #= filtered_subset_dds_DE
                                dds_res , #= subset_results
                                save_dir , #= subset_dir
                                padj_cutoff = 0.05,
                                annotation , #= annotation
                                name_file = name ){
  # Convert to table and add rownames to column with name "ensembl_id"
  results_table <- as_tibble(dds_res, rownames = "ensembl_id")
  # add annotation
  results_table <- left_join(results_table,annotation, by = c("ensembl_id" = "Gene.stable.ID"))
  # generate average counts per million metric from raw count data 
  raw_counts <- counts(dds_DE, normalized = F)
  cpms_table <- enframe(rowMeans(cpm(raw_counts)))
  colnames(cpms_table) <- c("ensembl_id", "avg_raw_cpm")
  # add avg cpms to table
  results_table <- left_join(results_table,cpms_table, by = c("ensembl_id" = "ensembl_id"))
  
  # Remove all NA from gene.name because this can give error when fetching entrez_id
  results_table <- results_table[!is.na(results_table$Gene.name),]
  # Add tag if gene is a DE gene
  results_table <- mutate(results_table, DE_gene = case_when(padj <= padj_cutoff ~ TRUE, T ~ FALSE))
  # add entrez id from GSEA
  results_table[(results_table$Gene.name == ""),"Gene.name"] <- "Not_found"
  # Get entrez ids from database
  entrez_ids <- results_table$Gene.name
  names(entrez_ids) <- mget(results_table$Gene.name,revmap(org.Mm.egSYMBOL),ifnotfound = NA)
  results_table[["Entrez_GSEA"]] <- names(entrez_ids) 
  # check for NA columns
  colSums(is.na(results_table))
  # Store all unfiltered data
  all_data <- results_table
  
  # Remove empty string gene names ("")
  results_table <- results_table[(results_table$Gene.name != "Not_found"),]
  # Remove duplicate gene names
  results_table <- results_table[!duplicated(results_table$Gene.name), ]
  # Remove non protein coding genes
  results_table <- results_table[(results_table$Gene.type == "protein_coding"),]
  # Remove all NA from padj 
  results_table <- results_table[!is.na(results_table$padj),]
  colSums(is.na(results_table))
  # Remove genes with different entrez id's
  dim(results_table)
  results_table <- results_table[(results_table$Entrez.id == results_table$Entrez_GSEA),]
  dim(results_table)
  # Store all filtered data
  filtered_all_data <- na.omit(results_table)
  colSums(is.na(filtered_all_data))
  
  
  dir.create(save_dir, showWarnings = FALSE)
  write.csv (all_data, file = paste0(save_dir,"/",name_file,".csv"), row.names = F)
  write.csv (filtered_all_data, file = paste0(save_dir,"/",name_file,"_filtered.csv"), row.names = F)
  write_xlsx(all_data, path = paste0(save_dir,"/",name_file,".xlsx"))
  write_xlsx(filtered_all_data, path = paste0(save_dir,"/",name_file,"_filtered.xlsx"))
}

# VolcanoPlots
get_paper_volcano_plot <- function(result_table = results,
                                   important_genes,
                                   title="test",
                                   crop_value_cutoff = 1e-10,
                                   padj_cutoff = 0.05,
                                   logFC_cutoff = 1,
                                   gene_num = 10,
                                   colors = c("#00FFFF", "gray50", "#ff1493"),
                                   point_size = 2.5,
                                   text_title_size = 15,
                                   x_scale_breaks = c(2,4,6,8,10),
                                   x_axis_limits = c(-12, 12) ){
  # Adding upregulated or downregulated tag
  fulldata <- result_table
  result_table <- mutate(result_table,Expression = case_when(log2FoldChange >= logFC_cutoff & padj <= padj_cutoff ~ "Up-regulated",
                                                             log2FoldChange <= -(logFC_cutoff) & padj <= padj_cutoff ~ "Down-regulated",
                                                             T ~ "Unchanged"))
  
  result_table <- mutate(result_table,Newpvalue = case_when(pvalue > crop_value_cutoff  ~ pvalue,
                                                            T ~ crop_value_cutoff))
  
  result_table <- mutate(result_table,View = case_when(pvalue > crop_value_cutoff  ~ "Outside",
                                                       T ~ "Inside"))
  
  fulldata <- mutate(fulldata,Expression = case_when(log2FoldChange >= logFC_cutoff & padj <= padj_cutoff ~ colors[3],
                                                     log2FoldChange <= -(logFC_cutoff) & padj <= padj_cutoff ~ colors[1],
                                                     T ~ colors[2]))
  
  fulldata <- mutate(fulldata,Newpvalue = case_when(pvalue > crop_value_cutoff  ~ pvalue,
                                                    T ~ crop_value_cutoff))
  
  fulldata <- mutate(fulldata,View = case_when(pvalue > crop_value_cutoff  ~ "Outside",
                                               T ~ "Inside"))
  
  # Filter on padj before doing this then order and take top X
  significant_genes <- result_table[result_table$padj < padj_cutoff,]
  significant_genes_ordered <- significant_genes[order(significant_genes$log2FoldChange,decreasing = T),]
  up <- head(significant_genes_ordered,round(gene_num/2,0))
  down <- tail(significant_genes_ordered,round(gene_num/2,0))
  top_genes <- rbind(up,down)
  top_genes_names <- top_genes$Gene.name
  
  # Get important genes
  found_genes <- c()
  for (gene in important_genes){
    if (gene %in% fulldata$Gene.name){
      found_genes <- append(found_genes,gene)
    }
  }
  data <- result_table[result_table$Gene.name == "nothing",]
  for (gene in found_genes){
    gene_info <- fulldata[fulldata$Gene.name == gene,]
    data <- rbind(data,gene_info)
  }
  top_genes_names <- data$Gene.name
  #print(top_genes_names)
  expr <- data$Expression
  #print(expr)
  
  volcano_plot <- ggplot(result_table, aes(log2FoldChange, -log(Newpvalue,10))) + # -log10 conversion  
    geom_point(aes(color = Expression, shape = View) , size = point_size) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values=c(1, 16), labels = c("Outside","Inside"))+ guides(shape = "none") +
    geom_label_repel(data = data,
                     mapping = aes(log2FoldChange, -log(Newpvalue,10), label = top_genes_names ,),
                     max.overlaps = 50,
                     size = point_size + 4) +
    #geom_vline(xintercept = c(-(logFC_cutoff),logFC_cutoff), linetype = "dashed", color = "black", size = 0.5) +
    ggtitle(title) +
    xlab(expression("log"[2]*"FC")) + 
    ylab(expression("-log"[10]*"p-value")) +
    scale_x_continuous(breaks = x_scale_breaks) +
    #coord_fixed(xlim=x_axis_limits) + 
    theme_classic() +
    theme(axis.title.x = element_text(face = "bold", color = "black",size = text_title_size-2,hjust = 0.5),
          axis.title.y = element_text(face = "bold", color = "black",size = text_title_size-2,hjust = 0.5), 
          axis.text.x = element_text(color = "black", size = text_title_size-8),
          axis.text.y = element_text(color = "black", size = text_title_size-8),
          plot.title = element_text(face = "bold", color = "black",size = text_title_size-2,hjust = 0.5),
          legend.text = element_text(size = text_title_size-8),
          legend.title = element_text(size = text_title_size-7, face = "bold"),
          legend.position = c(0.2,0.92)
          #legend.position=legend_pos,
          #legend.box.background = element_rect(colour = "black")
    )
  print(volcano_plot)
}

# GSEA PATHWAY PLOT
GSEA_plot <- function (fgsea_results,
                       NES_cutoff = 1.5,
                       npathw = 20,
                       color_down = "#00adef",
                       color_up = "#ff1493",
                       title = "No title",
                       reverse_order = T,
                       y_text_size = 7 ) {
  n <- npathw/2
  fgsea_results <- mutate(fgsea_results, short_name = str_split_fixed(fgsea_results$pathway,"_",2)[,2])
  fgsea_results$short_name <- str_replace_all(fgsea_results$short_name,"_"," ")
  print(fgsea_results$short_name)
  fgsea_results$short_name <- str_replace_all(fgsea_results$short_name," UP","")
  fgsea_results$short_name <- str_replace_all(fgsea_results$short_name," DN","")
  # Add table tag
  fgsea_results <- mutate(fgsea_results, Regulated = case_when(NES > 0 ~ "Up-regulated",
                                                               NES < 0 ~ "Down-regulated",
                                                               T ~ "Unchanged"))
  # order rows
  if (reverse_order == T){
    # Fetch top 20
    up <- fgsea_results[fgsea_results$Regulated == "Up-regulated",]
    up <- up[order(up$NES, decreasing = F),]
    
    down <- fgsea_results[fgsea_results$Regulated == "Down-regulated",]
    down <- down[order(down$NES, decreasing = F),]
    
    top20 <- rbind(tail(up,n),head(down,n))
  } else {
    up <- fgsea_results[fgsea_results$Regulated == "Up-regulated",]
    up <- up[order(up$NES, decreasing = F),]
    
    down <- fgsea_results[fgsea_results$Regulated == "Down-regulated",]
    down <- down[order(down$NES, decreasing = F),]
    
    top20 <- rbind(head(down,n),tail(up,n))
  }
  
  # Check for short names already present and remove them for now
  top20 <- top20[duplicated(top20$short_name) == F,]
  # lock in factor level order
  top20$short_name <- factor(top20$short_name, levels = top20$short_name)
  
  p <- ggplot(top20,aes(short_name, NES)) +
    geom_bar(stat= "identity", aes(fill = Regulated))+
    scale_fill_manual(values=c(color_down, color_up)) +
    coord_flip() +
    #scale_y_continuous(limits = c(-2, 2)) +
    geom_hline(yintercept = 0) +
    labs(x = "", y = "NES")+
    theme(text = element_text(family = "Calibri",color = "black"),
          #axis.text.y = element_text(size = y_text_size), 
          #plot.title = element_text(hjust = 1),
          axis.ticks.y = element_blank(),
          axis.line.x = element_line(),
          panel.background = element_rect(fill = "white"),
          #axis.text = element_text(color = "black"),
          axis.title = element_text(face = "plain", color = "black",size = 28,hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 25),
          axis.text.x = element_text(color = "black", size = 25),
          plot.title = element_text(face = "bold", color = "black",size = 26,hjust = 0.5), #0.5
          legend.position = "None",
          #legend.position = "top",
          #legend.text = element_text(size = 23),
          #legend.title = element_text(size = 26, face = "bold")
    ) +
    ggtitle(title)
  print(p)
}

# PCA plots
get_PCA <- function(vsd_obj,title){
  pcaData <- plotPCA(vsd_obj, intgroup=c("shape", "color"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  pca <- ggplot(pcaData, aes(PC1, PC2, color=shape, shape=color)) +
    geom_point(size=5) +
    scale_shape_manual(values=c(17, 15),name = NULL) +
    scale_color_manual(values=c("Blue","Red"),name = NULL) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    #ggtitle(title) +
    coord_fixed() +
    theme_classic() + 
    theme(#legend.position=c(.89, .15),
          legend.title = element_blank(),
          legend.key.size = unit(0.5, "cm"),
          legend.spacing.y = unit(0.0001, 'cm'),
          text = element_text(size = 20))
  print(pca) + stat_ellipse(type = "norm", linetype = 2, level = 0.8)
  
}