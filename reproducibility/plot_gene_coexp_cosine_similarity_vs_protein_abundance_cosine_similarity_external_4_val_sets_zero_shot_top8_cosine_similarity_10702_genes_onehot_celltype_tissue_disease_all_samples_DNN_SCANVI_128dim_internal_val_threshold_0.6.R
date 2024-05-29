#-----------------------------------------------------------------------------------------------------------------------
#External validation:

library(dplyr)
library(scales)
library(ggplot2)

#-----------------------------------------------------------------------------------------------------------------------
#internal validation:
#On server, Run function:
#-----------------------------------------------------------------------------------------------------------------------
#GSE128639:
softmax_df_gene_coexp_cosine_similarity <- read.csv('GSE128639_softmax_df_gene_coexp_cosine_similarity_20240504.csv', stringsAsFactors = F, row.names = 1, check.names = F)
softmax_df_protein_coexp_cosine_similarity <- read.csv('GSE128639_softmax_df_protein_exp_cosine_similarity_20240504.csv', stringsAsFactors = F, row.names = 1, check.names = F)
dataset_id = 'GSE128639'
cor.test(as.numeric(unlist(softmax_df_gene_coexp_cosine_similarity)), as.numeric(unlist(softmax_df_protein_coexp_cosine_similarity))) #corr: 0.6200932 

pdf(paste0('plot_gene_coexp_cosine_similarity_vs_protein_abundance_cosine_similarity_', dataset_id, '_20240502.pdf'))

plot(as.numeric(unlist(softmax_df_gene_coexp_cosine_similarity)), as.numeric(unlist(softmax_df_protein_coexp_cosine_similarity)), xlab = 'Pairwise similarity of protein abundance', ylab = 'Pairwise similarity of gene coexpression')
dev.off()

#-----------------------------------------------------------------------------------------------------------------------
#GSM5025059:
softmax_df_gene_coexp_cosine_similarity <- read.csv('GSM5025059_softmax_df_gene_coexp_cosine_similarity_20240504.csv', stringsAsFactors = F, row.names = 1, check.names = F)
softmax_df_protein_coexp_cosine_similarity <- read.csv('GSM5025059_softmax_df_protein_exp_cosine_similarity_20240504.csv', stringsAsFactors = F, row.names = 1, check.names = F)
dataset_id = 'GSM5025059'
cor.test(as.numeric(unlist(softmax_df_gene_coexp_cosine_similarity)), as.numeric(unlist(softmax_df_protein_coexp_cosine_similarity))) #corr: 0.7018749 

pdf(paste0('plot_gene_coexp_cosine_similarity_vs_protein_abundance_cosine_similarity_', dataset_id, '_20240502.pdf'))

plot(as.numeric(unlist(softmax_df_gene_coexp_cosine_similarity)), as.numeric(unlist(softmax_df_protein_coexp_cosine_similarity)), xlab = 'Pairwise similarity of protein abundance', ylab = 'Pairwise similarity of gene coexpression')
dev.off()
#-----------------------------------------------------------------------------------------------------------------------
#GSM5025052:
softmax_df_gene_coexp_cosine_similarity <- read.csv('GSM5025052_softmax_df_gene_coexp_cosine_similarity_20240504.csv', stringsAsFactors = F, row.names = 1, check.names = F)
softmax_df_protein_coexp_cosine_similarity <- read.csv('GSM5025052_softmax_df_protein_exp_cosine_similarity_20240504.csv', stringsAsFactors = F, row.names = 1, check.names = F)
dataset_id = 'GSM5025052'
cor.test(as.numeric(unlist(softmax_df_gene_coexp_cosine_similarity)), as.numeric(unlist(softmax_df_protein_coexp_cosine_similarity))) #corr: 0.749338 

pdf(paste0('plot_gene_coexp_cosine_similarity_vs_protein_abundance_cosine_similarity_', dataset_id, '_20240502.pdf'))

plot(as.numeric(unlist(softmax_df_gene_coexp_cosine_similarity)), as.numeric(unlist(softmax_df_protein_coexp_cosine_similarity)), xlab = 'Pairwise similarity of protein abundance', ylab = 'Pairwise similarity of gene coexpression')
dev.off()

#-----------------------------------------------------------------------------------------------------------------------
#GSM5093918:
softmax_df_gene_coexp_cosine_similarity <- read.csv('GSM5093918_softmax_df_gene_coexp_cosine_similarity_20240504.csv', stringsAsFactors = F, row.names = 1, check.names = F)
softmax_df_protein_coexp_cosine_similarity <- read.csv('GSM5093918_softmax_df_protein_exp_cosine_similarity_20240504.csv', stringsAsFactors = F, row.names = 1, check.names = F)
dataset_id = 'GSM5093918'
cor.test(as.numeric(unlist(softmax_df_gene_coexp_cosine_similarity)), as.numeric(unlist(softmax_df_protein_coexp_cosine_similarity))) #corr: 0.5218956 

pdf(paste0('plot_gene_coexp_cosine_similarity_vs_protein_abundance_cosine_similarity_', dataset_id, '_20240502.pdf'))

plot(as.numeric(unlist(softmax_df_gene_coexp_cosine_similarity)), as.numeric(unlist(softmax_df_protein_coexp_cosine_similarity)), xlab = 'Pairwise similarity of protein abundance', ylab = 'Pairwise similarity of gene coexpression')
dev.off()








