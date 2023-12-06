library(dplyr)
library(Seurat)
library(ggplot2)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

#-------------------------------------------------------------
#Use ADT data for cell-cell interaction analysis:
#
#Import data:
load('CT_HCC_scRNA_data.RData')

DefaultAssay(data.combined) = "RNA"
data.combined <- subset(data.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mito < 30)
features = rownames(data.combined)[rowSums(data.combined@assays$RNA@counts > 0) >= 5]
data.combined <- subset(data.combined, features = features)
data.combined <- NormalizeData(data.combined, normalization.method = "LogNormalize", scale.factor = 10000)

#data.combined <- FindVariableFeatures(data.combined)
data.combined <- ScaleData(data.combined)
#data.combined <- RunPCA(data.combined, verbose = FALSE)
#data.combined <- FindNeighbors(data.combined, dims = 1:50)
#data.combined <- FindClusters(data.combined, resolution = cluster_resolution)
#data.combined <- RunUMAP(data.combined, dims = 1:50) 

#Impute surface proteins as predicted:
#Seen proteins:
seen_all_protein_list <- read.csv('protein_gene_names_union_289_DNNs_from_combined_6_training_sets_DNN_onehot_celltype_tissue_disease_SCANVI_128dim_internal_val_threshold_0.6_20230115.csv', stringsAsFactors = F)
unseen_all_protein_list <- read.csv('max_softmax_Zi_Zt_zero_shot_top_8_10702_gene_cell_features_from_HCC_gdT_control_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_unseen_proteins_20230115.csv', stringsAsFactors = F)$X
unseen_all_protein_list = setdiff(unseen_all_protein_list, seen_all_protein_list$gene_name) #2554 unseen proteins (excluding overlap from 289 seen proteins)
unseen_all_protein_list = setdiff(unseen_all_protein_list, seen_all_protein_list$consistent_protein_name)

#Import seen proteins' predicted values:  
seen_file_pred_control <- NULL
seen_file_pred_tumor <- NULL

for (i in seen_all_protein_list$consistent_protein_name) {
  
  tmp_control = read.csv(paste0( i, '_cell_features_protein_specific_onehot_DNN_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_training_combined_6_sets_y_pred_seen_HCC_gdT_control_universal_features_SCANVI_128dim_impute_HVG_zeros_20230115.csv'), stringsAsFactors = F, row.names = 1)
  tmp_tumor = read.csv(paste0( i, '_cell_features_protein_specific_onehot_DNN_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_training_combined_6_sets_y_pred_seen_HCC_gdT_tumor_universal_features_SCANVI_128dim_impute_HVG_zeros_20230115.csv'), stringsAsFactors = F, row.names = 1)

  colnames(tmp_control)[1] = i
  colnames(tmp_tumor)[1] = i
  
  if(length(seen_file_pred_control) == 0){
    seen_file_pred_control = tmp_control
  }else{
    seen_file_pred_control = cbind(seen_file_pred_control, tmp_control)
  }
  
  if(length(seen_file_pred_tumor) == 0){
    seen_file_pred_tumor = tmp_tumor
  }else{
    seen_file_pred_tumor = cbind(seen_file_pred_tumor, tmp_tumor)
  }
  
}

seen_file_pred_tumor = t(seen_file_pred_tumor[, seen_all_protein_list$consistent_protein_name]) #289 16709
seen_file_pred_control = t(seen_file_pred_control[, seen_all_protein_list$consistent_protein_name]) #289 22211
seen_file_pred = cbind(seen_file_pred_tumor, seen_file_pred_control)[, colnames(data.combined)]
##matching_protein_gene_table <- read.csv('protein_gene_names_union_289_DNNs_from_combined_6_training_sets_DNN_onehot_celltype_tissue_disease_SCANVI_128dim_internal_val_threshold_0.8_20230115.csv', stringsAsFactors = F)
##rownames(matching_protein_gene_table) = matching_protein_gene_table$consistent_protein_name
##rownames(seen_file_pred) = matching_protein_gene_table[rownames(seen_file_pred), 'gene_name']

#Import unseen proteins' prdicted values:  
unseen_file_pred_control <- NULL
unseen_file_pred_tumor <- NULL

for (i in unseen_all_protein_list) {
  
  tmp_control = read.csv(paste0('zero_shot_top8_Zi_Zt_cosine_similarity_10702_gene_cell_features_from_HCC_gdT_control_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_pred_vs_truth_y_unseen_', i, '_20230115.csv'), stringsAsFactors = F, row.names = 1)
  tmp_tumor = read.csv(paste0('zero_shot_top8_Zi_Zt_cosine_similarity_10702_gene_cell_features_from_HCC_gdT_tumor_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_pred_vs_truth_y_unseen_', i, '_20230115.csv'), stringsAsFactors = F, row.names = 1)

  colnames(tmp_control)[1] = i
  colnames(tmp_tumor)[1] = i
  
  if(length(unseen_file_pred_control) == 0){
    unseen_file_pred_control = tmp_control
  }else{
    unseen_file_pred_control = cbind(unseen_file_pred_control, tmp_control)
  }
  
  if(length(unseen_file_pred_tumor) == 0){
    unseen_file_pred_tumor = tmp_tumor
  }else{
    unseen_file_pred_tumor = cbind(unseen_file_pred_tumor, tmp_tumor)
  }
  
}

unseen_file_pred_tumor = t(unseen_file_pred_tumor[, unseen_all_protein_list]) #289 16709
unseen_file_pred_control = t(unseen_file_pred_control[, unseen_all_protein_list]) #289 22211
unseen_file_pred = cbind(unseen_file_pred_tumor, unseen_file_pred_control)[, colnames(data.combined)]

#Select seen and unseen proteins with high confidence:
#
#filter proteins by internal validation performance:
seen_internal_val_performance <- read.csv('cor_per_pro_internal_val_6_combined_training_sets_cell_features_protein_specific_DNN_onehot_celltype_tissue_disease_SCANVI_128dim_64_32_16_0.0001_seen_proteins_20230115.csv', stringsAsFactors = F, row.names = 1)
use_seen_proteins <- rownames(seen_internal_val_performance)[seen_internal_val_performance$pearson > 0.6]

unseen_max_cosine_similarity_control <- read.csv('max_softmax_Zi_Zt_zero_shot_top_8_10702_gene_cell_features_from_HCC_gdT_control_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_unseen_proteins_20230115.csv', stringsAsFactors = F, row.names = 1)
unseen_max_cosine_similarity_tumor <- read.csv('max_softmax_Zi_Zt_zero_shot_top_8_10702_gene_cell_features_from_HCC_gdT_tumor_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_unseen_proteins_20230115.csv', stringsAsFactors = F, row.names = 1)

use_unseen_proteins <- intersect(rownames(unseen_max_cosine_similarity_control)[unseen_max_cosine_similarity_control$max_inferred_coef > 0.85], rownames(unseen_max_cosine_similarity_tumor)[unseen_max_cosine_similarity_tumor$max_inferred_coef > 0.85])

seen_file_pred = seen_file_pred[intersect(rownames(seen_file_pred), use_seen_proteins), ]
unseen_file_pred = unseen_file_pred[intersect(rownames(unseen_file_pred), use_unseen_proteins), ]

#
#For CellChat analysis, all the protein names should be changed to gene names:
#Check duplicated gene names for multiple proteins:
new_row_names <- rownames(seen_file_pred)

for (n in 1:nrow(seen_file_pred)) {
  
  if(length(intersect(rownames(match_training_protein_gene_name), rownames(seen_file_pred)[n])) > 0){
    if(match_training_protein_gene_name[rownames(seen_file_pred)[n] , 'gene_name'] != ' '){
      new_row_names[n] = match_training_protein_gene_name[rownames(seen_file_pred)[n] , 'gene_name']
    }
  }
  
}

rownames(seen_file_pred)[which(duplicated(new_row_names))]
seen_file_pred = seen_file_pred[setdiff(rownames(seen_file_pred), c("CD45", "CD45RO", "CD45RB", "HLA-ABC", "Folate")), ] #Manually check
#
for (n in 1:nrow(seen_file_pred)) {
  
  if(length(intersect(rownames(match_training_protein_gene_name), rownames(seen_file_pred)[n])) > 0){
    if(match_training_protein_gene_name[rownames(seen_file_pred)[n] , 'gene_name'] != ' '){
      rownames(seen_file_pred)[n] = match_training_protein_gene_name[rownames(seen_file_pred)[n] , 'gene_name']
    }
  }
  
}

seen_unseeen_file_pred <- rbind(seen_file_pred, unseen_file_pred) #394 proteins

#Replace proteins with negative predicted values with zeros, because CellChat cannot input negative data:
seen_unseeen_file_pred[seen_unseeen_file_pred < 0] = 0

data.combined[["ADT"]] <- CreateAssayObject(counts = seen_unseeen_file_pred)
data.combined@assays$ADT@data = seen_unseeen_file_pred
data.combined <- ScaleData(data.combined, assay = "ADT")

#-------------------------------------------------------------
#Cell-cell interaction analysis on imputed ADT:
DefaultAssay(data.combined) = "ADT"

data.combined.Tumor <- data.combined[, data.combined@meta.data$stim == 'Tumor']
data.combined.Control <- data.combined[, data.combined@meta.data$stim == 'Control']

cellchat.Tumor.ADT <- createCellChat(object = data.combined.Tumor, meta = data.combined.Tumor@meta.data, group.by = "predicted.id", assay = 'ADT')
cellchat.Control.ADT <- createCellChat(object = data.combined.Control, meta = data.combined.Control@meta.data, group.by = "predicted.id", assay = 'ADT')

for (i in c('cellchat.Tumor.ADT', 'cellchat.Control.ADT')) {
  
  cellchat <- get(i)
  CellChatDB <- CellChatDB.human 
  CellChatDB.use <- CellChatDB
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  #future::plan("multiprocess", workers = 4) # do parallel
  cellchat <- identifyOverExpressedGenes(cellchat, thresh.pc = 0.1, thresh.fc = 0.25, thresh.p = 0.05)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  assign(i, cellchat)
  
}

object.list <- list(Control = cellchat.Control.ADT, Tumor = cellchat.Tumor.ADT)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#
#Plot numbers / weights of interactions or interaction strength among different cell populations:
celltype_tumor_group = as.data.frame(table(data.combined@meta.data$predicted.id[data.combined@meta.data$stim == 'Tumor']))
celltype_tumor_group <- as.character(celltype_tumor_group$Var1[celltype_tumor_group$Freq >=3])
celltype_control_group <- as.data.frame(table(data.combined@meta.data$predicted.id[data.combined@meta.data$stim == 'Control']))
celltype_control_group <- as.character(celltype_control_group$Var1[celltype_control_group$Freq >=3])

celltype = intersect(celltype_tumor_group, celltype_control_group)

for(c in celltype){
  
pdf(paste0(c, '_cell_cell_interaction_counts_use_predicted_ADT_cellchat_HCC_gdT_20230115.pdf'))
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  tmp <- object.list[[i]]@net$count
  tmp[setdiff(rownames(tmp), c), ] = 0
  if(sum(tmp) != 0){
  netVisual_circle(tmp, weight.scale = T, label.edge= T, edge.weight.max = weight.max[2], edge.width.max = 4, edge.label.cex = 0.4, vertex.label.cex = 0.4, title.name = paste0("Number of interactions - ", names(object.list)[i]))
  }
}
dev.off()

#
pdf(paste0(c, 'cell_cell_interaction_weights_use_predicted_ADT_cellchat_HCC_gdT_20230115.pdf'))
weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  tmp <- object.list[[i]]@net$weight
  tmp[setdiff(rownames(tmp), c), ] = 0
  if(sum(tmp) != 0){
  netVisual_circle(tmp, weight.scale = T, label.edge= T, edge.weight.max = weight.max[2], edge.width.max = 4, edge.label.cex = 0.4, vertex.label.cex = 0.4, title.name = paste0("Weight of interactions - ", names(object.list)[i]))
  }
}
dev.off()

}


























