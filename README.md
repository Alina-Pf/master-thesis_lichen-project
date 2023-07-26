# master-thesis_lichen-project
 A collection of scripts I created for my masters thesis.
# 
++ script: TE_blast_hits_on_kmer_PCA_automated.R
+ Consists of two separate loops:
  1. Plots blast hits on PCA based on "short" kmer (min length = 10'000)
      + also exports PCA tables
  2. Plots GC content on PCA based on "long" kmer (min length = 20'000) (with short kmers, GC content is "averaged"/no differenced)
#
++ script: compare_all_lists_with_all_lists
 + Input: two lists listing names of files to be compared
#
++ cript: count_KO_identifiers.R
+ creates a table with:
+   number of total proteins
+   number of SP proteins (without tmd)
+   number of KOs detected
+   number of enriched/multiplied KOs
