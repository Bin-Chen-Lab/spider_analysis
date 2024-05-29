#Step 1: Gather PPI for all proteins in the combined 6 training set.
#Step 2: Replace current gene co-expression with PPI.
#Step 3: Obtain internal validation performance for unseen proteins.
#2021-8-22
#Use PPI data
library(dplyr)
library(STRINGdb)

#On PC:
#-----------------------------------------------------------------------------------------------------------------------------------------
#changeable parameters:
#GSE128639:
protein_gene_feature <- read.csv('all_25_protein_gene_name_GSE128639_BM_20230115.csv', stringsAsFactors = FALSE, check.names = F)
protein_gene_feature_emsemble_members <- read.csv('gene_names_selected_166_unseen_proteins_from_combined_6_training_sets_DNN_internal_val_20230115.csv', stringsAsFactors = FALSE, check.names = F)

dataset_id = 'GSE128639'
#-----------------------------------------------------------------------------------------------------------------------------------------
#changeable parameters:
#GSM5025059:
protein_gene_feature <- read.csv('all_13_protein_gene_name_GSM5025059_20230115.csv', stringsAsFactors = FALSE, check.names = F)
protein_gene_feature_emsemble_members <- read.csv('gene_names_selected_166_unseen_proteins_from_combined_6_training_sets_DNN_internal_val_20230115.csv', stringsAsFactors = FALSE, check.names = F)

dataset_id = 'GSM5025059'
#-----------------------------------------------------------------------------------------------------------------------------------------
#changeable parameters:
#GSM5025052:
protein_gene_feature <- read.csv('all_13_protein_gene_name_GSM5025052_20230115.csv', stringsAsFactors = FALSE, check.names = F)
protein_gene_feature_emsemble_members <- read.csv('gene_names_selected_166_unseen_proteins_from_combined_6_training_sets_DNN_internal_val_20230115.csv', stringsAsFactors = FALSE, check.names = F)

dataset_id = 'GSM5025052'
#-----------------------------------------------------------------------------------------------------------------------------------------
#changeable parameters:
#GSM5093918:
protein_gene_feature <- read.csv('all_39_protein_gene_name_GSM5093918_20230115.csv', stringsAsFactors = FALSE, check.names = F)
protein_gene_feature_emsemble_members <- read.csv('gene_names_selected_166_unseen_proteins_from_combined_6_training_sets_DNN_internal_val_20230115.csv', stringsAsFactors = FALSE, check.names = F)

dataset_id = 'GSM5093918'
#-----------------------------------------------------------------------------------------------------------------------------------------
gene_symbols <- unique(c(protein_gene_feature$gene_name, protein_gene_feature_emsemble_members$gene_name)) #166 genes


landmark_genes <- data.frame(gene_name = gene_symbols)
string_db <- STRINGdb$new( version="11.5", species=9606,
                           score_threshold=0, input_directory="")
#'https://stringdb-static.org/download/protein.info.v11.5/9606.protein.info.v11.5.txt.gz'
#'https://stringdb-static.org/download/protein.links.v11.5/9606.protein.links.v11.5.txt.gz'
#string_db$get_interactions(c('9606.ENSP00000000412', '9606.ENSP00000001008'))
#STRINGdb$help("map")
data_mapped <- string_db$map(landmark_genes, my_data_frame_id_col_names = 'gene_name', removeUnmappedRows = TRUE)
data_mapped = data_mapped[-which(duplicated(data_mapped$gene_name)), ]
write.csv(data_mapped, paste0('mapping_gene_name_STRINGid_external_val_set_', dataset_id, '.csv'))

#On server:
#-----------------------------------------------------------------------------------------------------------------------------------------
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
data_mapped <- read.csv(paste0('mapping_gene_name_STRINGid_external_val_set_', dataset_id, '.csv'), stringsAsFactors = F, header = T, check.names = F)
data_links <- read.table('9606.protein.links.v11.5.txt', stringsAsFactors = F, header = T)
data_links = data_links[data_links$combined_score >= 400, ]
length(unique(data_links$protein1))
length(unique(data_links$protein2))

#Construct the PPI matrix:
PPI_matrix <- matrix(0, nrow = length(unique(data_links$protein1)), ncol = length(unique(data_links$protein1)))
rownames(PPI_matrix) = unique(data_links$protein1)
colnames(PPI_matrix) = unique(data_links$protein1)

for (n in 1:nrow(data_links)) {
  
  PPI_matrix[data_links$protein1[n], data_links$protein2[n]] = data_links$combined_score[n]

}

PPI_matrix = PPI_matrix[data_mapped$STRING_id, ]
rownames(PPI_matrix) = data_mapped$gene_name

#delete columns contaning only zeros:
PPI_matrix = PPI_matrix[, colSums(PPI_matrix) > 0]

#Remove ensemble members that doesn't have PPI:
use_ensemble_members = protein_gene_feature_emsemble_members$gene_name
use_ensemble_members = use_ensemble_members[use_ensemble_members != "IGHG1"]
use_ensemble_members = use_ensemble_members[use_ensemble_members != "IGHM"]
use_ensemble_members = use_ensemble_members[use_ensemble_members != "TRGV9"]
use_ensemble_members = use_ensemble_members[use_ensemble_members != "IGKC"]

test_unseen_proteins = unique(protein_gene_feature$gene_name)
if(dataset_id == "GSM5093918"){
  test_unseen_proteins = test_unseen_proteins[test_unseen_proteins != "TRG"]
  test_unseen_proteins = test_unseen_proteins[test_unseen_proteins != "TRAV1-2"]
  test_unseen_proteins = test_unseen_proteins[test_unseen_proteins != "TRAV24"]
}

write.csv(PPI_matrix[test_unseen_proteins, ], paste0(dataset_id, '_external_val_set_protein_human_PPI_combined_score_400_threshold_20240502.csv'), quote=FALSE)
write.csv(PPI_matrix[use_ensemble_members, ], paste0(dataset_id, '_external_val_set_protein_human_PPI_combined_score_400_threshold_corresponding_ensemble_members_20240502.csv'), quote=FALSE)











