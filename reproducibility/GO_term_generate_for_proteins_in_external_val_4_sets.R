#Step 1: Gather GO terms for all proteins in the combined 6 training set.
#Step 2: Replace current gene co-expression with GO terms.
#Step 3: Obtain internal validation performance for unseen proteins.

library(dplyr)
library(ontologySimilarity)

data(gene_GO_terms) # Gene Ontology annotation of genes

#-----------------------------------------------------------------------------------------------------------------
#changeable parameters:
#GSE128639:
protein_gene_feature <- read.csv('all_25_protein_gene_name_GSE128639_BM_20230115.csv', stringsAsFactors = FALSE, check.names = F)
protein_gene_feature_emsemble_members <- read.csv('gene_names_selected_166_unseen_proteins_from_combined_6_training_sets_DNN_internal_val_20230115.csv', stringsAsFactors = FALSE, check.names = F)

dataset_id = 'GSE128639'
#-----------------------------------------------------------------------------------------------------------------
#changeable parameters:
#healthy pancreas GSM5025059:
protein_gene_feature <- read.csv('all_13_protein_gene_name_GSM5025059_20230115.csv', stringsAsFactors = FALSE, check.names = F)
protein_gene_feature_emsemble_members <- read.csv('gene_names_selected_166_unseen_proteins_from_combined_6_training_sets_DNN_internal_val_20230115.csv', stringsAsFactors = FALSE, check.names = F)

dataset_id = 'GSM5025059'
#-----------------------------------------------------------------------------------------------------------------
#changeable parameters:
#pancreas_pancreatitis_GSM5025052:
protein_gene_feature <- read.csv('all_13_protein_gene_name_GSM5025052_20230115.csv', stringsAsFactors = FALSE, check.names = F)
protein_gene_feature_emsemble_members <- read.csv('gene_names_selected_166_unseen_proteins_from_combined_6_training_sets_DNN_internal_val_20230115.csv', stringsAsFactors = FALSE, check.names = F)

dataset_id = 'GSM5025052'
#-----------------------------------------------------------------------------------------------------------------
#changeable parameters:
#COVID BALF GSM5093918:
protein_gene_feature <- read.csv('all_39_protein_gene_name_GSM5093918_20230115.csv', stringsAsFactors = FALSE, check.names = F)
protein_gene_feature_emsemble_members <- read.csv('gene_names_selected_166_unseen_proteins_from_combined_6_training_sets_DNN_internal_val_20230115.csv', stringsAsFactors = FALSE, check.names = F)

dataset_id = 'GSM5093918'
#-----------------------------------------------------------------------------------------------------------------
gene_symbols <- unique(c(protein_gene_feature$gene_name, protein_gene_feature_emsemble_members$gene_name)) #166 genes
landmark_genes_GO_terms <- gene_GO_terms[gene_symbols]

# gather all go terms for the lm genes
all_landmark_genes_go_terms_vector <- c()
for (i in (1:length(gene_symbols) )){
  go_selected <- landmark_genes_GO_terms[[i]]
  all_landmark_genes_go_terms_vector <- c(all_landmark_genes_go_terms_vector, go_selected)
}
all_unique_landmark_genes_go_terms <- sort(unique(all_landmark_genes_go_terms_vector)) #1551 go terms

# get go terms that are attributed to by more than 0 landmark genes
common_go_terms <- c()
for (i in (1:length(all_unique_landmark_genes_go_terms))){
  lm_go_term<-all_unique_landmark_genes_go_terms[[i]]
  count <- 0
  for (j in (1:length(gene_symbols))){
    lm_gene_go_terms <- landmark_genes_GO_terms[[j]]
    if (any(lm_gene_go_terms==lm_go_term)){
      count <- count+1
    }
    #if (count > 3){
    if (count > 0){
      common_go_terms <- c(common_go_terms, lm_go_term)
      break
    }
  }
}

# descriptors binary occurrence vector
gene_go_fingerprint <- matrix(0,length(gene_symbols),length(common_go_terms))
for (i in (1:length(common_go_terms))){
  for (j in (1:length(gene_symbols))){
    if (length(which(landmark_genes_GO_terms[[j]] == common_go_terms[[i]])) == 1){
      gene_go_fingerprint[j,i] <- 1
    }
  }
}

rownames(gene_go_fingerprint) <- gene_symbols
colnames(gene_go_fingerprint) <- common_go_terms
write.csv(gene_go_fingerprint[unique(protein_gene_feature$gene_name), ], paste0('combined_6_training_sets_protein_gene_GO_fingerprints_at_least_1_apperance_', dataset_id, '_20240428.csv'), quote=FALSE)
write.csv(gene_go_fingerprint[protein_gene_feature_emsemble_members$gene_name, ], paste0('combined_6_training_sets_protein_gene_GO_fingerprints_at_least_1_apperance_', dataset_id, '_corresponding_ensemble_members_20240428.csv'), quote=FALSE)




