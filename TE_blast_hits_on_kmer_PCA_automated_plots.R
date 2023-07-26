#########
# One run for all plots
#########

library(ggplot2)
library(colorspace)
library(tibble)

##### list with names of samples
ls <- read.table("~/data/dir_TE/list_of_data_sets", quote="\"", comment.char="")
ls <- as.list(ls)

##### named vector for colors
#myPalette <- choose_palette()
#myPalette(6) # prints hexcodes for the colors

class_cols <- c(" Bacteria " = "#C87A8A",
                " Algae " = "#45A271", # A29048
                " Fungi " = "#A29048", # 45A271
                " Bacteria  Algae " = "#2A9EB5",
                " Fungi  Algae " = "#A782C3")
class_cols


##### loop to produce all plots [TE his]
for (i in 1:length(ls[["V1"]])) {
  samp_name <- ls[["V1"]][i]
  setwd("~/data/dir_TE/dir_assemblies")
  data_kmer <- read.delim(paste0("wicKmer_output_", samp_name, "_size_500_10000-10000000"))
  
  setwd("~/data/dir_TE/dir_log_blast_strong_hit_eval")
  data_TE <- read.delim(paste0("log_blast_strong_hit_ctg_eval_", samp_name, "-100"))
  
  as.factor(data_TE$class)
  data_c <- merge(data_TE, data_kmer, by = "ctg", all.y = TRUE)
  data_c <- data.frame(data_c, row.names=data_c[,1])
  
  pca_df <- data_c[6:261]
  pca <- prcomp(pca_df, scale. = TRUE)
  
  pca_export <- as.data.frame(pca$x)
  pca_export <- rownames_to_column(pca_export, "ctg")
  pca_export <- pca_export[1:4]
  
  
  PCi_new <- data.frame(pca$x, class = data_c$class, hits = data_c$hits)
  PC_vals <- summary(pca)$importance[2,]
  PC1_val <- 100*PC_vals[1]
  PC2_val <- 100*PC_vals[2]
  
  the_ggplot <- ggplot(PCi_new, aes(x=PC1, y=PC2)) +
    geom_point(alpha=0.1, size = 0.7) +
    #geom_point(aes(x=PC1, y=PC2, col=class), size=log(PCi_new$hits)) +
    geom_point(aes(x=PC1, y=PC2, col=class), size=log(3*PCi_new$hits)) +
    scale_colour_manual(values = class_cols) +
    ggtitle(paste0("PCA of kmer (k=4) analysis of ", samp_name,"\nColored by TE blast hits, hit lenght = 100")) +
    xlab(paste0("PC1 (",PC1_val,"%)")) + ylab(paste0("PC2 (", PC2_val,"%)")) +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 23),
          plot.title = element_text(size = 27),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 23))
  
  # "small" kmers (from min size 10'000) -> GC content evenly everywhere
  the_ggplot2 <- ggplot(PCi_new, aes(x=PC1, y=PC2, color = data_kmer$GC)) +
    geom_point(size = 0.8) +
    ggtitle(paste0("PCA of kmer (k=4) analysis of ", samp_name,"\nColored by GC content")) +
    xlab(paste0("PC1 (",PC1_val,"%)")) + ylab(paste0("PC2 (", PC2_val,"%)")) +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 23),
          plot.title = element_text(size = 27),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 23)) +
    labs(color='GC')
  
  # export what is needed
  setwd("~/data/dir_TE/dir_R_plot_loop_output")
  #write.table(pca_export, sep= "\t", file=paste0("PCA_table_", samp_name), row.names = FALSE)
  #ggsave(the_ggplot, filename = paste0("TE_blast_hits_on_kmer_PCA_", samp_name, "-100.png"), width = 34, height = 22, units = "cm")
  #ggsave(the_ggplot2, filename = paste0("kmer_PCA_", samp_name, "_GC.png"), width = 34, height = 22, units = "cm")
}

################################################################################

### GC kmer plots, "long" wicKmer contigs (min length 20'000)

ls <- read.table("~/data/dir_lichensamples_collected_2022/dir_wicKmer_output/list_of_names_kmer_data_sets", quote="\"", comment.char="")
ls <- as.list(ls)

for (i in 1:length(ls[["V1"]])) {
  samp_name <- ls[["V1"]][i]
  setwd("~/data/dir_lichensamples_collected_2022/dir_wicKmer_output")
  data_1 <- read.delim(paste0("wicKmer_output_", samp_name, "_50000-10000000"))
  data_2 <- read.delim(paste0("wicKmer_output_", samp_name, "_20000-50000"))
  data_kmer <- rbind(data_1, data_2)
  
  pca_df <- data_kmer[3:259]
  pca <- prcomp(pca_df, scale. = TRUE)
  
  PCi_new <- data.frame(pca$x)
  PC_vals <- summary(pca)$importance[2,]
  PC1_val <- 100*PC_vals[1]
  PC2_val <- 100*PC_vals[2]
  
  the_ggplot <- ggplot(PCi_new, aes(x=PC1, y=PC2, color = data_kmer$GC)) +
    geom_point(size = 0.8) +
    ggtitle(paste0("PCA of kmer (k=4) analysis of ", samp_name,"\nColored by GC content")) +
    xlab(paste0("PC1 (",PC1_val,"%)")) + ylab(paste0("PC2 (", PC2_val,"%)")) +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 23),
          plot.title = element_text(size = 27),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 23)) +
    labs(color='GC')
  
  setwd("~/data/dir_TE/dir_R_plot_loop_output/dir_PCA_GC")
  #write.table(pca_export, sep= "\t", file=paste0("PCA_table_", samp_name), row.names = FALSE)
  ggsave(the_ggplot, filename = paste0("kmer_PCA_", samp_name, "_GC.png"), width = 34, height = 22, units = "cm")
}

