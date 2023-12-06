library(dplyr)
library(Seurat)
library(ggplot2)
#-------------------------------------------------------------
#Import RNA data:
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
unseen_all_protein_list <- read.csv('zero_shot_learning_protein_specific_onehot_celltype_tissue_disease_DNN_ensemble/max_softmax_Zi_Zt_zero_shot_top_8_10702_gene_cell_features_from_HCC_gdT_control_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_unseen_proteins_20230115.csv', stringsAsFactors = F)$X
unseen_all_protein_list = setdiff(unseen_all_protein_list, seen_all_protein_list$gene_name) #2554 unseen proteins (excluding overlap from 289 seen proteins)
unseen_all_protein_list = setdiff(unseen_all_protein_list, seen_all_protein_list$consistent_protein_name)

#Import seen proteins' prdicted values:  
seen_file_pred_control <- NULL
seen_file_pred_tumor <- NULL

for (i in seen_all_protein_list$consistent_protein_name) {
    
    tmp_control = read.csv(paste0(i, '_cell_features_protein_specific_onehot_DNN_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_training_combined_6_sets_y_pred_seen_HCC_gdT_control_universal_features_SCANVI_128dim_impute_HVG_zeros_20230115.csv'), stringsAsFactors = F, row.names = 1)
    tmp_tumor = read.csv(paste0(i, '_cell_features_protein_specific_onehot_DNN_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_training_combined_6_sets_y_pred_seen_HCC_gdT_tumor_universal_features_SCANVI_128dim_impute_HVG_zeros_20230115.csv'), stringsAsFactors = F, row.names = 1)

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

#Import unseen proteins' prdicted values:  
unseen_file_pred_control <- NULL
unseen_file_pred_tumor <- NULL

for (i in unseen_all_protein_list) {
  
  tmp_control = read.csv(paste0('zero_shot_top8_Zi_Zt_cosine_similarity_10702_gene_cell_features_from_HCC_gdT_control_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_pred_vs_truth_y_unseen_', i, '_20230115.csv'), stringsAsFactors = F, row.names = 1)
  tmp_tumor = read.csv(paste0('zero_shot_learning_protein_specific_onehot_celltype_tissue_disease_DNN_ensemble/zero_shot_top8_Zi_Zt_cosine_similarity_10702_gene_cell_features_from_HCC_gdT_tumor_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_pred_vs_truth_y_unseen_', i, '_20230115.csv'), stringsAsFactors = F, row.names = 1)

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

seen_unseeen_file_pred <- rbind(seen_file_pred, unseen_file_pred)

data.combined[["ADT"]] <- CreateAssayObject(counts = seen_unseeen_file_pred)
data.combined@assays$ADT@data = seen_unseeen_file_pred
data.combined <- ScaleData(data.combined, assay = "ADT")
  
#############################################################################################
#Analysis 1:
#Find DEGs between HCC control group and tumor group for each cell type, respectively:
celltype_tumor_group = as.data.frame(table(data.combined@meta.data$predicted.id[data.combined@meta.data$stim == 'Tumor']))
celltype_tumor_group <- as.character(celltype_tumor_group$Var1[celltype_tumor_group$Freq >=3])
celltype_control_group <- as.data.frame(table(data.combined@meta.data$predicted.id[data.combined@meta.data$stim == 'Control']))
celltype_control_group <- as.character(celltype_control_group$Var1[celltype_control_group$Freq >=3])

celltype = intersect(celltype_tumor_group, celltype_control_group)

for (c in celltype) {
  tmp <- data.combined[, colnames(data.combined)[data.combined@meta.data$predicted.id == c]]
  rna_markers <- FindMarkers(tmp, ident.1 = colnames(tmp)[tmp@meta.data$stim == 'Tumor'], ident.2 = colnames(tmp)[tmp@meta.data$stim == 'Control'], assay = "RNA")
  adt_markers <- FindMarkers(tmp, ident.1 = colnames(tmp)[tmp@meta.data$stim == 'Tumor'], ident.2 = colnames(tmp)[tmp@meta.data$stim == 'Control'], assay = "ADT")
  assign(paste0(c, '_DEGs'), rna_markers)
  assign(paste0(c, '_DEPs'), adt_markers)
  
}

rm(tmp)
rm(tmp_tumor)
rm(tmp_control)
save.image('case_study_HCC_gdT_DEGs_DEPs.RData')
#-------------------------------------------------------------------------------------------------------
#Use heatmaps to show DEPs and DEGs, respectively. Using p values < 0.05 to filter DEPs and DEGs:
#filter proteins by internal validation performance:
seen_internal_val_performance <- read.csv('cor_per_pro_internal_val_6_combined_training_sets_cell_features_protein_specific_DNN_onehot_celltype_tissue_disease_SCANVI_128dim_64_32_16_0.0001_seen_proteins_20230115.csv', stringsAsFactors = F, row.names = 1)
use_seen_proteins <- rownames(seen_internal_val_performance)[seen_internal_val_performance$pearson > 0.6]

unseen_max_cosine_similarity_control <- read.csv('zero_shot_learning_protein_specific_onehot_celltype_tissue_disease_DNN_ensemble/max_softmax_Zi_Zt_zero_shot_top_8_10702_gene_cell_features_from_HCC_gdT_control_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_unseen_proteins_20230115.csv', stringsAsFactors = F, row.names = 1)
unseen_max_cosine_similarity_tumor <- read.csv('zero_shot_learning_protein_specific_onehot_celltype_tissue_disease_DNN_ensemble/max_softmax_Zi_Zt_zero_shot_top_8_10702_gene_cell_features_from_HCC_gdT_tumor_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_unseen_proteins_20230115.csv', stringsAsFactors = F, row.names = 1)

use_unseen_proteins <- intersect(rownames(unseen_max_cosine_similarity_control)[unseen_max_cosine_similarity_control$max_inferred_coef > 0.88], rownames(unseen_max_cosine_similarity_tumor)[unseen_max_cosine_similarity_tumor$max_inferred_coef > 0.88])


for (c in celltype) {
  
  DEPs <- filter(get(paste0(c, '_DEPs')), p_val_adj < 0.05)
  DEGs <- filter(get(paste0(c, '_DEGs')), p_val_adj < 0.05)
  
  
  use_protein = intersect(c(use_seen_proteins, use_unseen_proteins), rownames(DEPs))
  use_gene = rownames(DEGs)
  
  #if(length(use_gene) > 30){
  #  use_gene = rownames(get(paste0(c, '_DEGs')))[order(abs(get(paste0(c, '_DEGs'))$avg_log2FC), decreasing = T)[1:30]]
  #}
  
  #if(length(use_protein) > 30){
   # use_protein = rownames(get(paste0(c, '_DEPs')))[order(abs(get(paste0(c, '_DEPs'))$avg_log2FC), decreasing = T)[1:30]]
  #}
  
  tmp <- data.combined[, colnames(data.combined)[data.combined@meta.data$predicted.id == c]]
  
  if(length(use_gene) > 0){
    
  pdf(paste0('heatmap_DEGs_', c, '.pdf'), height = 8, width = 8)

  p_DEGs = DoHeatmap(object = tmp, assay = 'RNA', features = use_gene, cells = colnames(tmp), group.by = "stim", disp.min = 0, disp.max = 1) + 
    #scale_fill_gradientn(colors = c("blue", "white", "red"))+
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 8, name = "RdBu")))+
    theme(text = element_text(size = 10))
    
  print(p_DEGs)
  dev.off()
  }
  
  if(length(use_protein) > 0){
  
  pdf(paste0('heatmap_DEPs_', c, '.pdf'), height = 8, width = 8)

  p_DEPs = DoHeatmap(object = tmp, assay = 'ADT', features = use_protein, cells = colnames(tmp), group.by = "stim", disp.min = 0, disp.max = 1) + 
      #scale_fill_gradientn(colors = c("blue", "white", "red"))+
      scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 8, name = "RdBu")))+
      theme(text = element_text(size = 10))
  print(p_DEPs)
  dev.off()
  }
  
}

#-------------------------------------------------------------------------------------------------------
#Compare DEGs and DEPs to find difference. Use Featureplot to show the proteins and genes that were not overlapped:
match_training_protein_gene_name = read.csv('protein_gene_names_union_289_DNNs_from_combined_6_training_sets_DNN_onehot_celltype_tissue_disease_SCANVI_128dim_internal_val_threshold_0.6_20230115.csv', stringsAsFactors = F)
rownames(match_training_protein_gene_name) = match_training_protein_gene_name$consistent_protein_name

for (c in celltype) {
  
  DEPs <- filter(get(paste0(c, '_DEPs')), p_val_adj < 0.05)
  DEGs <- filter(get(paste0(c, '_DEGs')), p_val_adj < 0.05)
  
  use_protein = intersect(c(use_seen_proteins, use_unseen_proteins), rownames(DEPs))
  use_gene = rownames(DEGs)
  diff_DEGs_DEPs <- setdiff(setdiff(use_protein, rownames(match_training_protein_gene_name)), use_gene)
  overlap_DEGs_DEPs <- intersect(setdiff(use_protein, rownames(match_training_protein_gene_name)), use_gene)
  
  for (n in 1:length(use_protein)) {
    
    if(length(intersect(rownames(match_training_protein_gene_name), use_protein[n])) > 0){
      if(length(intersect(match_training_protein_gene_name[use_protein[n] , 'gene_name'], use_protein)) == 0){
        if(match_training_protein_gene_name[use_protein[n] , 'gene_name'] != ' '){
          if(length(setdiff(match_training_protein_gene_name[use_protein[n] , 'gene_name'], use_gene)) > 0){
            diff_DEGs_DEPs = c(diff_DEGs_DEPs, use_protein[n])
          }else{
            overlap_DEGs_DEPs = c(overlap_DEGs_DEPs, use_protein[n])
    }
        }else{
          diff_DEGs_DEPs = c(diff_DEGs_DEPs, use_protein[n])
      }
      }else{
        diff_DEGs_DEPs = c(diff_DEGs_DEPs, use_protein[n])
    }
    }
  }

  assign(paste0(c, '_diff_DEGs_DEPs'), unique(diff_DEGs_DEPs))
  assign(paste0(c, '_overlap_DEGs_DEPs'), unique(overlap_DEGs_DEPs))
  
}

#-------------------------------------------------------------------------------------------------------
#Use Featureplot to show the proteins and genes that were not overlapped:
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
      if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
       hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
color_list <- ggplotColours(n=18)
set.seed(1)
color_list = sample(color_list)

pdf('visualize_case_study_HCC_gdT_by_celltype_20230115.pdf', width = 12, height = 10)
p = DimPlot(data.combined, reduction = "umap",label=F,group.by='predicted.id',
            cols = color_list)
print(p)
dev.off()

pdf('visualize_case_study_HCC_gdT_by_patient_20230115.pdf', width = 10.5, height = 10)
p = DimPlot(data.combined, reduction = "umap",label=F,group.by='orig.ident')
print(p)
dev.off()

pdf('visualize_case_study_HCC_gdT_by_disease_20230115.pdf', width = 10.5, height = 10)
p = DimPlot(data.combined, reduction = "umap",label=F,group.by='stim')
print(p)
dev.off()


#-------------------------------------------------------------------------------------------------------------------
#Plot FeaturePlot of RNA and protein expressions:
for (c in celltype) {
  
  all_protein_list <- rownames(get(paste0(c, '_DEPs')))[get(paste0(c, '_DEPs'))$avg_log2FC < (-0.8) | get(paste0(c, '_DEPs'))$avg_log2FC > 0.8]
  all_protein_list = intersect(all_protein_list, get(paste0(c, '_diff_DEGs_DEPs')))
  
  if(length(all_protein_list) > 0){
  for (p in all_protein_list) {
    
    #plot protein:
    pdf(paste0('Predicted_', c, '_', p, '.pdf'), height = 5, width = 6)
    DefaultAssay(data.combined) <- "ADT"
    p1 = FeaturePlot(data.combined, features = p, min.cutoff = "q05", max.cutoff = "q95")
    print(p1)
    dev.off()
    
    #plot gene:
    if(length(intersect(rownames(match_training_protein_gene_name), p)) > 0){
        if(match_training_protein_gene_name[p , 'gene_name'] != ' '){
            gene = match_training_protein_gene_name[p , 'gene_name']
        }else{
          gene = p
        }
    }else{
      gene = p
    }
    pdf(paste0('RNA_', c, '_', p, '.pdf'), height = 5, width = 6)
    DefaultAssay(data.combined) <- "RNA"
    p1 = FeaturePlot(data.combined, features = gene, min.cutoff = "q05", max.cutoff = "q95")
    print(p1)
    dev.off()
    
  }
  }
}

#
for (c in celltype) {
  
  all_protein_list <- intersect(rownames(get(paste0(c, '_DEPs')))[order(get(paste0(c, '_DEPs'))$avg_log2FC, decreasing = T)], c(use_seen_proteins, use_unseen_proteins))[1:10]
  all_protein_list <- c(all_protein_list, intersect(rownames(get(paste0(c, '_DEPs')))[order(get(paste0(c, '_DEPs'))$avg_log2FC, decreasing = F)], c(use_seen_proteins, use_unseen_proteins)))[1:10]

  if(length(all_protein_list) > 0){
    for (p in all_protein_list) {
      
      #plot protein:
      pdf(paste0('Predicted_', c, '_', p, '_top10FC.pdf'), height = 5, width = 6)
      DefaultAssay(data.combined) <- "ADT"
      p1 = FeaturePlot(data.combined, features = p, min.cutoff = "q05", max.cutoff = "q95")
      print(p1)
      dev.off()
      
      #plot gene:
      if(length(intersect(rownames(match_training_protein_gene_name), p)) > 0){
        if(match_training_protein_gene_name[p , 'gene_name'] != ' '){
          gene = match_training_protein_gene_name[p , 'gene_name']
        }else{
          gene = p
        }
      }else{
        gene = p
      }
      
      if(length(intersect(gene, rownames(data.combined@assays$RNA))) > 0){
      pdf(paste0('RNA_', c, '_', p, '_top10FC.pdf'), height = 5, width = 6)
      DefaultAssay(data.combined) <- "RNA"
      p1 = FeaturePlot(data.combined, features = gene, min.cutoff = "q05", max.cutoff = "q95")
      print(p1)
      dev.off()
      }
      
    }
  }
}


#############################################################################################

#############################################################################################
#Analysis 2:
#Find DEGs and DEPs for each cell type:
load('case_study_HCC_gdT_DEGs_DEPs.RData')
DefaultAssay(data.combined) = "RNA"

celltype = as.data.frame(table(data.combined@meta.data$predicted.id))

#filter proteins by internal validation performance:
seen_internal_val_performance <- read.csv('cor_per_pro_internal_val_6_combined_training_sets_cell_features_protein_specific_DNN_onehot_celltype_tissue_disease_SCANVI_128dim_64_32_16_0.0001_seen_proteins_20230115.csv', stringsAsFactors = F, row.names = 1)
use_seen_proteins <- rownames(seen_internal_val_performance)[seen_internal_val_performance$pearson > 0.6]

unseen_max_cosine_similarity_control <- read.csv('max_softmax_Zi_Zt_zero_shot_top_8_10702_gene_cell_features_from_HCC_gdT_control_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_unseen_proteins_20230115.csv', stringsAsFactors = F, row.names = 1)
unseen_max_cosine_similarity_tumor <- read.csv('max_softmax_Zi_Zt_zero_shot_top_8_10702_gene_cell_features_from_HCC_gdT_tumor_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_unseen_proteins_20230115.csv', stringsAsFactors = F, row.names = 1)

use_unseen_proteins <- intersect(rownames(unseen_max_cosine_similarity_control)[unseen_max_cosine_similarity_control$max_inferred_coef > 0.85], rownames(unseen_max_cosine_similarity_tumor)[unseen_max_cosine_similarity_tumor$max_inferred_coef > 0.85])

match_training_protein_gene_name = read.csv('protein_gene_names_union_289_DNNs_from_combined_6_training_sets_DNN_onehot_celltype_tissue_disease_SCANVI_128dim_internal_val_threshold_0.6_20230115.csv', stringsAsFactors = F)
rownames(match_training_protein_gene_name) = match_training_protein_gene_name$consistent_protein_name


tmp <- data.combined[, colnames(data.combined)]
  
#
#Find DEGs:
celltype_DE <- as.character(as.data.frame(table(data.combined@active.ident))$Var1)
use_celltype_DE = celltype_DE[-c(which(celltype_DE == 'Erythroid.cells'), which(celltype_DE == 'Hep.5'))] #Removing cell types with < 3 cells under either condition. These cells need to use FindAllMarkers() to identify DEGs / DEPs.

for (i in use_celltype_DE) {
  tmp2 <- FindConservedMarkers(tmp, ident.1 = i, grouping.var = "stim", assay = 'RNA', verbose = FALSE)
  tmp2$cluster = i
  assign(paste0('celltype_DEGs_', i), tmp2)
}

#
#Find DEPs:
tmp <- data.combined[, colnames(data.combined)]

use_protein = intersect(c(use_seen_proteins, use_unseen_proteins), rownames(tmp@assays$ADT))
use_protein_2 = use_protein

index <- NULL

for (n in 1:length(use_protein_2)) {
  
  if(length(intersect(rownames(match_training_protein_gene_name), use_protein_2[n])) > 0){
    if(length(intersect(match_training_protein_gene_name[use_protein_2[n] , 'gene_name'], use_protein_2)) > 1){
      if(match_training_protein_gene_name[use_protein_2[n] , 'gene_name'] != ' '){
        if(match_training_protein_gene_name[use_protein_2[n] , 'gene_name'] != use_protein_2[n]){
          index <- c(index, which(use_protein_2 ==  match_training_protein_gene_name[use_protein_2[n] , 'gene_name']))
        }else{
          index <- c(index, which(use_protein_2 ==  match_training_protein_gene_name[use_protein_2[n] , 'gene_name']))[2]
        }
      }
    }
  }
  
  if(length(intersect(rownames(match_training_protein_gene_name), use_protein_2[n])) > 0){
    if(length(intersect(match_training_protein_gene_name[use_protein_2[n] , 'gene_name'], use_protein_2)) == 1){
      if(match_training_protein_gene_name[use_protein_2[n] , 'gene_name'] != ' '){
        if(match_training_protein_gene_name[use_protein_2[n] , 'gene_name'] != use_protein_2[n]){
          index <- c(index, which(use_protein_2 ==  match_training_protein_gene_name[use_protein_2[n] , 'gene_name']))
        }
      }
    }
  }
}

if(length(index) > 0){
  use_protein_2 = use_protein_2[-index]
}

use_protein_2 = use_protein_2[!is.na(use_protein_2)]

DefaultAssay(tmp) = "ADT"
tmp <- subset(tmp, features = use_protein_2)

for (i in use_celltype_DE) {
  tmp2 <- FindConservedMarkers(tmp, ident.1 = i, grouping.var = "stim", assay = 'ADT', verbose = FALSE)
  tmp2$cluster = i
  assign(paste0('celltype_DEPs_', i), tmp2)
}


tmp <- data.combined[, colnames(data.combined)]

celltype_DEGs_all <- NULL

for (i in use_celltype_DE) {
  
  tmp2 <- get(paste0('celltype_DEGs_', i))
  tmp2$gene = rownames(tmp2)
  
  if(length(celltype_DEGs_all) == 0){
    celltype_DEGs_all = tmp2
  }else{
    celltype_DEGs_all = rbind(celltype_DEGs_all, tmp2)
  }
  
}

celltype_DEGs_all <-  celltype_DEGs_all[celltype_DEGs_all$minimump_p_val < 0.05, ]
celltype_DEGs_all %>% group_by(cluster) %>% slice_max(n = 2, order_by = Control_avg_log2FC)
celltype_DEGs_all %>% group_by(cluster) %>% top_n(n = 10, wt = Control_avg_log2FC) -> top10

#
#DEG heatmap:
pdf(paste0('heatmap_conserved_DEGs_combine_tumor_and_control_20230824.pdf'), height = 8, width = 8)

p_DEGs = DoHeatmap(object = tmp, assay = 'RNA', features = top10$gene, disp.min = 0, disp.max = 1) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 8, name = "RdYlGn")))+
  theme(text = element_text(size = 3))

print(p_DEGs)

dev.off()

#
#DEP heatmap:
DefaultAssay(tmp) = "ADT"
tmp <- subset(tmp, features = use_protein_2)

celltype_DEPs_all <- NULL

for (i in use_celltype_DE) {
  
  tmp2 <- get(paste0('celltype_DEPs_', i))
  tmp2$protein = rownames(tmp2)
  
  if(length(celltype_DEPs_all) == 0){
    celltype_DEPs_all = tmp2
  }else{
    celltype_DEPs_all = rbind(celltype_DEPs_all, tmp2)
  }
  
}


celltype_DEPs_all <-  celltype_DEPs_all[celltype_DEPs_all$minimump_p_val < 0.05, ]
celltype_DEPs_all %>% group_by(cluster) %>% slice_max(n = 2, order_by = Control_avg_log2FC)
celltype_DEPs_all %>% group_by(cluster) %>% top_n(n = 10, wt = Control_avg_log2FC) -> top10

pdf(paste0('heatmap_conserved_DEPs_combine_tumor_and_control_20230824.pdf'), height = 8, width = 8)

p_DEPs = DoHeatmap(object = tmp, assay = 'ADT', features = top10$protein, disp.min = 0, disp.max = 1) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 8, name = "RdYlGn")))+
  theme(text = element_text(size = 3))

print(p_DEPs)

dev.off()









