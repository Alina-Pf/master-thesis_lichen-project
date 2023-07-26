#########
# Count KO numbers and duplicated KO numbers for each sample!
#########

library(dplyr)


### prep outside of loop
df_temp=data.frame()
df_export=data.frame()

#ls <- read.table("~/data/dir_secretome/dir_genomes_to_analyze/dir_SP_removed_descr_lines/list_of_sample_names_without_Neurospora", quote="\"", comment.char="")
#ls <- read.table("~/data/dir_secretome/dir_genomes_to_analyze/list_of_sample_names_for_R", quote="\"", comment.char="")
ls <- read.table("~/data/dir_secretome/dir_sp_prot_for_alina/dir_KOassign_output/list_of_genomes_transmembrane_rm", quote="\"", comment.char="")
ls <- as.list(ls)


for (i in 1:length(ls[["V1"]])) {
  genome_name <- ls[["V1"]][i]
  
  ### all proteins
  setwd("~/data/dir_secretome/dir_genomes_to_analyze/dir_descr_lines_ALL_proteins")
  descr_lines_file <- read.delim(paste0(genome_name,"_CDS_proteins_descr_lines_reformatted"), header=FALSE)
  total_protein <- nrow(descr_lines_file)
  
  ### only SP proteins
  #setwd("~/data/dir_secretome/dir_KO_input_for_R")
  setwd("~/data/dir_secretome/dir_sp_prot_for_alina/dir_KOassign_output")
  df <- read.csv(paste0("user_ko_",genome_name,".csv"), header=FALSE)
  total_SP <- nrow(df)
  
  ### all detected KOs [KEGG Orthologs]
  df <- replace(df, df == "", NA)
  df$V2 <- as.factor(df$V2)
  count_KOs <- nlevels(df$V2)
  all_KOs <- levels(df$V2)
  
  ### duplicated KOs
  facLevel <- lapply(df, table)
  test <- sapply(facLevel, function(x) x[x != 1])
  new_df <- test[["V2"]]
  new_df2 <- as.data.frame(new_df)
  new_df2$Freq <- as.integer(new_df2$Freq)
  
  duplicated_KOs <- nrow(new_df2)
  KO_identifiers <- new_df2$Var1
  
  df_temp = cbind(genome_name, total_protein, total_SP, count_KOs, duplicated_KOs)
  
  ### output stuff
  df_export <- rbind(df_export, df_temp)
  #setwd("~/data/dir_secretome/dir_genomes_to_analyze/dir_R_output_genome_comparison")
  setwd("~/data/dir_secretome/dir_sp_prot_for_alina/dir_R_output")
  write.table(df_export, quote = FALSE, sep= "\t", row.names = FALSE, col.names = TRUE, file="genomes_comparison_overview_new")
  #write.table(new_df2, quote = FALSE, sep= "\t", row.names = FALSE, col.names = TRUE, file=paste0(genome_name,"_duplicated_KOss"))
  #write.table(KO_identifiers, quote = FALSE, sep= "\t", row.names = FALSE, col.names = FALSE, file=paste0(genome_name,"_duplicated_KO_identifiers"))
  write.table(all_KOs, quote = FALSE, sep= "\t", row.names = FALSE, col.names = FALSE, file=paste0(genome_name,"_all_KO_identifiers"))
}

#### single files manually #####################################################


### data handling
# description lines of all proteins
setwd("~/data/dir_secretome/dir_genomes_to_analyze/dir_descr_lines_ALL_proteins")
descr_lines_file <- read.delim("Neurospora_crassa_CDS_proteins_descr_lines_reformatted", header=FALSE)
total_protein <- nrow(descr_lines_file)

# KO data set, converted to .csv data type
setwd("~/data/dir_secretome/dir_KO_input_for_R")
df <- read.csv("user_ko_genome_R.csv", header=FALSE)

genome_name <- "loop_name"
head(df)

### count total SP proteins (=total rows)
total_SP <- nrow(df)


### count total KOs 
df <- replace(df, df == "", NA)

df$V2 <- as.factor(df$V2) 
count_KOs <- nlevels(df$V2)
count_KOs # correct, empty/NA is not counted as factor
all_KOs <- levels(df$V2)

### count duplicated KOS (facteor occurrence >1)
facLevel <- lapply(df, table)
test <- sapply(facLevel, function(x) x[x != 1])

new_df <- test[["V2"]]
head(new_df)

new_df2 <- as.data.frame(new_df)
head(new_df2)
new_df2$Freq <- as.integer(new_df2$Freq)

duplicated_KOs <- nrow(new_df2)
duplicated_KOs
KO_identifiers <- new_df2$Var1

duplicated_KOs_ex_separate <- new_df2$Var1
duplicated_KOs_ex_separate


### put total_SP, count_KOs, duplicated_KOs into ONE ROW
df_export = cbind(genome_name,total_protein,total_SP,count_KOs,duplicated_KOs)
df_export # this....works?????
# does rbind() work?
# does merge() work?


setwd("~/data/dir_secretome/dir_genomes_to_analyze/dir_R_output_genome_comparison")
#write.table(new_df2, quote = FALSE, sep= "\t", row.names = FALSE, col.names = TRUE, file="Neurospora_crassa_duplicated_KOs")
write.table(KO_identifiers, quote = FALSE, sep= "\t", row.names = FALSE, col.names = FALSE, file="genome_R_duplicated_KO_identifiers")
write.table(all_KOs, quote = FALSE, sep= "\t", row.names = FALSE, col.names = FALSE, file="genome_R_all_KO_identifiers")


