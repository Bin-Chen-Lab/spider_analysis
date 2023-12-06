library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)


#############################################################################################
#Analysis 1:
#Find DEGs / DEPs between different sites (tissues) for each cell type, respectively:
#primary_treatment_naïve vs. liver_metastases_treatment_naïve 
#primary_treatment_naïve vs. PBMC_treatment_naïve 
#liver_metastases_treatment_naïve vs. PBMC_treatment_naïve 
load('case_study_metastasis_GSE178318_DEGs_DEPs_20230826.RData')


table(filter(RNA_integrated@meta.data, study == 'PBMC_treatment_naïve'|study == 'PBMC_treated')$cell_type) #There is only 1 epithelial cell in PBMCs.

#change parameters:
group1 = 'primary_treatment_naïve'
group2 = 'liver_metastases_treatment_naïve'
#group1 = 'primary_treatment_naïve', group2 = 'PBMC_treatment_naïve'
#group1 = 'liver_metastases_treatment_naïve', group2 = 'PBMC_treatment_naïve'

celltype_group_1 = as.data.frame(table(RNA_integrated@meta.data$cell_type[RNA_integrated@meta.data$study == group1]))
celltype_group_1 <- as.character(celltype_group_1$Var1[celltype_group_1$Freq >=3])
celltype_group_2 <- as.data.frame(table(RNA_integrated@meta.data$cell_type[RNA_integrated@meta.data$study == group2]))
celltype_group_2 <- as.character(celltype_group_2$Var1[celltype_group_2$Freq >=3])

celltype = intersect(celltype_group_1, celltype_group_2)

for (c in celltype) {
  tmp <- RNA_integrated[, colnames(RNA_integrated)[RNA_integrated@meta.data$cell_type == c]]
  rna_markers <- FindMarkers(tmp, ident.1 = colnames(tmp)[tmp@meta.data$study == group1], ident.2 = colnames(tmp)[tmp@meta.data$study == group2], assay = "RNA")
  adt_markers <- FindMarkers(tmp, ident.1 = colnames(tmp)[tmp@meta.data$study == group1], ident.2 = colnames(tmp)[tmp@meta.data$study == group2], assay = "ADT")
  
  if(group1 == 'primary_treatment_naïve' & group2 == 'liver_metastases_treatment_naïve'){
    assign(paste0(c, '_DEGs_primary_TN_vs_LM_TN'), rna_markers)
    assign(paste0(c, '_DEPs_primary_TN_vs_LM_TN'), adt_markers)
  }
  
  if(group1 == 'primary_treatment_naïve' & group2 == 'PBMC_treatment_naïve'){
    assign(paste0(c, '_DEGs_primary_TN_vs_PBMC_TN'), rna_markers)
    assign(paste0(c, '_DEPs_primary_TN_vs_PBMC_TN'), adt_markers)
  }
  
  if(group1 == 'liver_metastases_treatment_naïve' & group2 == 'PBMC_treatment_naïve'){
    assign(paste0(c, '_DEGs_LM_TN_vs_PBMC_TN'), rna_markers)
    assign(paste0(c, '_DEPs_LM_TN_vs_PBMC_TN'), adt_markers)
  }
  
}

rm(tmp)

save.image('case_study_metastasis_GSE178318_DEGs_DEPs_20230826.RData')

#
#filter out highly-confident proteins:
load('case_study_metastasis_GSE178318_DEGs_DEPs_20230826.RData')

seen_internal_val_performance <- read.csv('cor_per_pro_internal_val_6_combined_training_sets_cell_features_protein_specific_DNN_onehot_celltype_tissue_disease_SCANVI_128dim_64_32_16_0.0001_seen_proteins_20230115.csv', stringsAsFactors = F, row.names = 1)
use_seen_proteins <- rownames(seen_internal_val_performance)[seen_internal_val_performance$pearson > 0.6]

for (d in c('primary_TN', 'LM_TN', 'PBMC_TN', 'primary_T', 'LM_T', 'PBMC_T')) {
  
  if(d == 'primary_TN'|d == 'LM_T'){
    tmp <- read.csv(paste0('max_softmax_Zi_Zt_zero_shot_top_8_10702_gene_cell_features_from_GSE178318_', d, '_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_unseen_proteins_20230115.csv'), stringsAsFactors = F, row.names = 1)
  }
  
  if(d == 'primary_T'){
    tmp <- read.csv(paste0('max_softmax_Zi_Zt_zero_shot_top_8_10699_gene_cell_features_from_GSE178318_', d, '_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_unseen_proteins_20230115.csv'), stringsAsFactors = F, row.names = 1)
  }
  
  if(d == 'LM_TN'){
    tmp <- read.csv(paste0('max_softmax_Zi_Zt_zero_shot_top_8_10698_gene_cell_features_from_GSE178318_', d, '_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_unseen_proteins_20230115.csv'), stringsAsFactors = F, row.names = 1)
  }
  
  if(d == 'PBMC_T'){
    tmp <- read.csv(paste0('max_softmax_Zi_Zt_zero_shot_top_8_10696_gene_cell_features_from_GSE178318_', d, '_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_unseen_proteins_20230115.csv'), stringsAsFactors = F, row.names = 1)
  }
  
  if(d == 'PBMC_TN'){
    tmp <- read.csv(paste0('max_softmax_Zi_Zt_zero_shot_top_8_10691_gene_cell_features_from_GSE178318_', d, '_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_unseen_proteins_20230115.csv'), stringsAsFactors = F, row.names = 1)
  }
  
  tmp = tmp[get(paste0('unseen_all_protein_list_', d)), ]
  tmp = as.data.frame(tmp, row.names = get(paste0('unseen_all_protein_list_', d)))
  colnames(tmp) = 'max_inferred_coef'
  assign(paste0('unseen_protein_max_similarity_', d), tmp)
  
}

#
#primary_treatment_naïve vs. liver_metastases_treatment_naïve: 
use_unseen_proteins <- intersect(rownames(unseen_protein_max_similarity_primary_TN)[unseen_protein_max_similarity_primary_TN$max_inferred_coef > 0.85], rownames(unseen_protein_max_similarity_LM_TN)[unseen_protein_max_similarity_LM_TN$max_inferred_coef > 0.85])

group1 = 'primary_treatment_naïve'
group2 = 'liver_metastases_treatment_naïve'
#group2 = 'PBMC_treatment_naïve'

celltype_group_1 = as.data.frame(table(RNA_integrated@meta.data$cell_type[RNA_integrated@meta.data$study == group1]))
celltype_group_1 <- as.character(celltype_group_1$Var1[celltype_group_1$Freq >=3])
celltype_group_2 <- as.data.frame(table(RNA_integrated@meta.data$cell_type[RNA_integrated@meta.data$study == group2]))
celltype_group_2 <- as.character(celltype_group_2$Var1[celltype_group_2$Freq >=3])

celltype = intersect(celltype_group_1, celltype_group_2)

#Compare DEGs and DEPs to find difference. Use Featureplot to show the proteins and genes that were not overlapped:
match_training_protein_gene_name = read.csv('protein_gene_names_union_289_DNNs_from_combined_6_training_sets_DNN_onehot_celltype_tissue_disease_SCANVI_128dim_internal_val_threshold_0.6_20230115.csv', stringsAsFactors = F)
rownames(match_training_protein_gene_name) = match_training_protein_gene_name$consistent_protein_name

DEPs_ratio <- NULL
DEGs_ratio <- NULL
select_DEGs_ratio <- NULL
select_DEPs_ratio <- NULL

for (c in celltype) {
  
  DEPs <- filter(get(paste0(c, '_DEPs_primary_TN_vs_LM_TN')), p_val_adj < 0.05)
  DEGs <- filter(get(paste0(c, '_DEGs_primary_TN_vs_LM_TN')), p_val_adj < 0.05)
  #DEPs <- filter(get(paste0(c, '_DEPs_primary_TN_vs_PBMC_TN')), p_val_adj < 0.05)
  #DEGs <- filter(get(paste0(c, '_DEGs_primary_TN_vs_PBMC_TN')), p_val_adj < 0.05)
  
  #-----------
  #Calculate the percentage of DEPs for each cell type, respectively:
  DEPs_ratio = c(DEPs_ratio, length(intersect(c(use_seen_proteins, use_unseen_proteins), rownames(DEPs))) / length(c(use_seen_proteins, use_unseen_proteins)))
  #Calculate the percentage of DEGs (use all genes) for each cell type, respectively:
  all_genes <- rownames(RNA_integrated@assays$RNA)
  DEGs_ratio = c(DEGs_ratio, length(rownames(DEGs)) / length(all_genes))
  #Calculate the percentage of DEGs (use proteins' corresponding genes) for each cell type, respectively:
  select_proteins <- NULL
  select_genes <- NULL
  for (p in c(use_seen_proteins)) {
    if(length(intersect(match_training_protein_gene_name[p, 'gene_name'], all_genes)) > 0 | length(intersect(p, all_genes)) > 0){
      select_proteins = c(select_proteins, p)
      if(length(intersect(match_training_protein_gene_name[p, 'gene_name'], all_genes)) > 0){
        select_genes = c(select_genes, match_training_protein_gene_name[p, 'gene_name'])
      }else if(length(intersect(p, all_genes)) > 0){
        select_genes = c(select_genes, p)
      }
    }
  }
  select_proteins = c(select_proteins, intersect(use_unseen_proteins, all_genes)) #287
  select_genes = c(select_genes, intersect(use_unseen_proteins, all_genes))
  select_DEGs_ratio = c(select_DEGs_ratio, length(intersect(select_genes, rownames(DEGs))) / length(select_genes))
  select_DEPs_ratio = c(select_DEPs_ratio, length(intersect(select_proteins, rownames(DEPs))) / length(select_proteins))
  #-----------
  
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
  
  assign(paste0(c, '_diff_DEGs_DEPs_primary_TN_vs_LM_TN'), unique(diff_DEGs_DEPs))
  assign(paste0(c, '_overlap_DEGs_DEPs_primary_TN_vs_LM_TN'), unique(overlap_DEGs_DEPs))
  #assign(paste0(c, '_diff_DEGs_DEPs_primary_TN_vs_PBMC_TN'), unique(diff_DEGs_DEPs))
  #assign(paste0(c, '_overlap_DEGs_DEPs_primary_TN_vs_PBMC_TN'), unique(overlap_DEGs_DEPs))
  
  
}

#-----------
#Calculate the percentage of DEPs & DEGs for each cell type, respectively:
d1 <- data.frame(cell_type = celltype, 'percentage_DEPs (100%)' = select_DEPs_ratio*100, stringsAsFactors = F)
d1 = d1[order(d1[, 2], decreasing = T), ]
colnames(d1) = c('cell_type', 'percentage_DEPs (100%)')
d2 <- data.frame(cell_type = celltype, 'percentage_DEGs (100%)' = select_DEGs_ratio*100, stringsAsFactors = F)
d2 = d2[order(d2[, 2], decreasing = T), ]
colnames(d2) = c('cell_type', 'percentage_DEGs (100%)')
d3 <- data.frame(cell_type = celltype, 'percentage_DEGs (100%)' = DEGs_ratio*100, stringsAsFactors = F)
d3 = d3[order(d3[, 2], decreasing = T), ]
colnames(d3) = c('cell_type', 'percentage_DEGs (100%)')
#------------

#-------------------------------------------------------------------------------------------------------
#Use Featureplot to show the proteins and genes that were not overlapped:
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
color_list <- ggplotColours(n=10)
set.seed(1)
color_list = sample(color_list)

pdf('visualize_case_study_GSE178318_by_celltype_20230903.pdf', width = 12, height = 10)
p = DimPlot(RNA_integrated, reduction = "umap",label=F,group.by='cell_type',
            cols = color_list, raster = F)
print(p)
dev.off()

#
#Plot FeaturePlot of RNA and protein expressions:
for (c in celltype) {
  
  all_protein_list = get(paste0(c, '_diff_DEGs_DEPs_primary_TN_vs_LM_TN'))
  
  if(length(all_protein_list) > 0){
    for (p in all_protein_list) {
      
      #plot protein:
      pdf(paste0('Predicted_', c, '_', p, '_primary_TN_vs_LM_TN_20230903.pdf'), height = 5, width = 6)
      DefaultAssay(RNA_integrated) <- "ADT"
      p1 = FeaturePlot(RNA_integrated, features = p, min.cutoff = "q05", max.cutoff = "q95", raster = F)
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
      pdf(paste0('RNA_', c, '_', p, '_primary_TN_vs_LM_TN_20230903.pdf'), height = 5, width = 6)
      DefaultAssay(RNA_integrated) <- "RNA"
      p1 = FeaturePlot(RNA_integrated, features = gene, min.cutoff = "q05", max.cutoff = "q95", raster = F)
      print(p1)
      dev.off()
      
    }
  }
}

#

#############################################################################################
#Analysis 2:
#Find DEGs / DEPs between different disease states (i.e., before / after chemotherapy) for each cell type, respectively:
load('case_study_metastasis_GSE178318_DEGs_DEPs_20230826.RData')

#change parameters:
#group1 = 'primary_treatment_naïve', group2 = 'primary_treated'
group1 = 'liver_metastases_treatment_naïve'
group2 = 'liver_metastases_treated'
#group1 = 'PBMC_treatment_naïve'
#group2 = 'PBMC_treated'


celltype_group_1 = as.data.frame(table(RNA_integrated@meta.data$cell_type[RNA_integrated@meta.data$study == group1]))
celltype_group_1 <- as.character(celltype_group_1$Var1[celltype_group_1$Freq >=3])
celltype_group_2 <- as.data.frame(table(RNA_integrated@meta.data$cell_type[RNA_integrated@meta.data$study == group2]))
celltype_group_2 <- as.character(celltype_group_2$Var1[celltype_group_2$Freq >=3])

celltype = intersect(celltype_group_1, celltype_group_2)

for (c in celltype) {
  tmp <- RNA_integrated[, colnames(RNA_integrated)[RNA_integrated@meta.data$cell_type == c]]
  rna_markers <- FindMarkers(tmp, ident.1 = colnames(tmp)[tmp@meta.data$study == group1], ident.2 = colnames(tmp)[tmp@meta.data$study == group2], assay = "RNA")
  adt_markers <- FindMarkers(tmp, ident.1 = colnames(tmp)[tmp@meta.data$study == group1], ident.2 = colnames(tmp)[tmp@meta.data$study == group2], assay = "ADT")
  
  if(group1 == 'primary_treatment_naïve' & group2 == 'primary_treated'){
    assign(paste0(c, '_DEGs_primary_TN_vs_primary_T'), rna_markers)
    assign(paste0(c, '_DEPs_primary_TN_vs_primary_T'), adt_markers)
  }
  
  if(group1 == 'liver_metastases_treatment_naïve' & group2 == 'liver_metastases_treated'){
    assign(paste0(c, '_DEGs_LM_TN_vs_LM_T'), rna_markers)
    assign(paste0(c, '_DEPs_LM_TN_vs_LM_T'), adt_markers)
  }
  
  if(group1 == 'PBMC_treatment_naïve' & group2 == 'PBMC_treated'){
    assign(paste0(c, '_DEGs_PBMC_TN_vs_PBMC_T'), rna_markers)
    assign(paste0(c, '_DEPs_PBMC_TN_vs_PBMC_T'), adt_markers)
  }
  
}

rm(tmp)

save.image('case_study_metastasis_GSE178318_DEGs_DEPs_20230826.RData')

#
seen_internal_val_performance <- read.csv('cor_per_pro_internal_val_6_combined_training_sets_cell_features_protein_specific_DNN_onehot_celltype_tissue_disease_SCANVI_128dim_64_32_16_0.0001_seen_proteins_20230115.csv', stringsAsFactors = F, row.names = 1)
use_seen_proteins <- rownames(seen_internal_val_performance)[seen_internal_val_performance$pearson > 0.6]

for (d in c('primary_TN', 'LM_TN', 'PBMC_TN', 'primary_T', 'LM_T', 'PBMC_T')) {
  
  if(d == 'primary_TN'|d == 'LM_T'){
    tmp <- read.csv(paste0('max_softmax_Zi_Zt_zero_shot_top_8_10702_gene_cell_features_from_GSE178318_', d, '_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_unseen_proteins_20230115.csv'), stringsAsFactors = F, row.names = 1)
  }
  
  if(d == 'primary_T'){
    tmp <- read.csv(paste0('max_softmax_Zi_Zt_zero_shot_top_8_10699_gene_cell_features_from_GSE178318_', d, '_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_unseen_proteins_20230115.csv'), stringsAsFactors = F, row.names = 1)
  }
  
  if(d == 'LM_TN'){
    tmp <- read.csv(paste0('max_softmax_Zi_Zt_zero_shot_top_8_10698_gene_cell_features_from_GSE178318_', d, '_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_unseen_proteins_20230115.csv'), stringsAsFactors = F, row.names = 1)
  }
  
  if(d == 'PBMC_T'){
    tmp <- read.csv(paste0('max_softmax_Zi_Zt_zero_shot_top_8_10696_gene_cell_features_from_GSE178318_', d, '_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_unseen_proteins_20230115.csv'), stringsAsFactors = F, row.names = 1)
  }
  
  if(d == 'PBMC_TN'){
    tmp <- read.csv(paste0('max_softmax_Zi_Zt_zero_shot_top_8_10691_gene_cell_features_from_GSE178318_', d, '_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_unseen_proteins_20230115.csv'), stringsAsFactors = F, row.names = 1)
  }
  
  tmp = tmp[get(paste0('unseen_all_protein_list_', d)), ]
  tmp = as.data.frame(tmp, row.names = get(paste0('unseen_all_protein_list_', d)))
  colnames(tmp) = 'max_inferred_coef'
  assign(paste0('unseen_protein_max_similarity_', d), tmp)
  
}

#
#primary_treatment_naïve vs. primary_treated: 
#use_unseen_proteins <- intersect(rownames(unseen_protein_max_similarity_primary_TN)[unseen_protein_max_similarity_primary_TN$max_inferred_coef > 0.85], rownames(unseen_protein_max_similarity_primary_T)[unseen_protein_max_similarity_primary_T$max_inferred_coef > 0.85])
use_unseen_proteins <- intersect(rownames(unseen_protein_max_similarity_LM_TN)[unseen_protein_max_similarity_LM_TN$max_inferred_coef > 0.85], rownames(unseen_protein_max_similarity_LM_T)[unseen_protein_max_similarity_LM_T$max_inferred_coef > 0.85])

#group1 = 'primary_treatment_naïve'
#group2 = 'primary_treated'
group1 = 'liver_metastases_treatment_naïve'
group2 = 'liver_metastases_treated'

celltype_group_1 = as.data.frame(table(RNA_integrated@meta.data$cell_type[RNA_integrated@meta.data$study == group1]))
celltype_group_1 <- as.character(celltype_group_1$Var1[celltype_group_1$Freq >=3])
celltype_group_2 <- as.data.frame(table(RNA_integrated@meta.data$cell_type[RNA_integrated@meta.data$study == group2]))
celltype_group_2 <- as.character(celltype_group_2$Var1[celltype_group_2$Freq >=3])

celltype = intersect(celltype_group_1, celltype_group_2)

#Compare DEGs and DEPs to find difference. Use Featureplot to show the proteins and genes that were not overlapped:
match_training_protein_gene_name = read.csv('protein_gene_names_union_289_DNNs_from_combined_6_training_sets_DNN_onehot_celltype_tissue_disease_SCANVI_128dim_internal_val_threshold_0.6_20230115.csv', stringsAsFactors = F)
rownames(match_training_protein_gene_name) = match_training_protein_gene_name$consistent_protein_name

for (c in celltype) {
  
  #DEPs <- filter(get(paste0(c, '_DEPs_primary_TN_vs_primary_T')), p_val_adj < 0.05)
  #DEGs <- filter(get(paste0(c, '_DEGs_primary_TN_vs_primary_T')), p_val_adj < 0.05)
  DEPs <- filter(get(paste0(c, '_DEPs_LM_TN_vs_LM_T')), p_val_adj < 0.05)
  DEGs <- filter(get(paste0(c, '_DEGs_LM_TN_vs_LM_T')), p_val_adj < 0.05)
  
  
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
  
  #assign(paste0(c, '_diff_DEGs_DEPs_primary_TN_vs_primary_T'), unique(diff_DEGs_DEPs))
  #assign(paste0(c, '_overlap_DEGs_DEPs_primary_TN_vs_primary_T'), unique(overlap_DEGs_DEPs))
  assign(paste0(c, '_diff_DEGs_DEPs_LM_TN_vs_LM_T'), unique(diff_DEGs_DEPs))
  assign(paste0(c, '_overlap_DEGs_DEPs_LM_TN_vs_LM_T'), unique(overlap_DEGs_DEPs))
  
  
}

#-------------------------------------------------------------------------------------------------------
#Use Featureplot to show the proteins and genes that were not overlapped:
#
#Plot FeaturePlot of RNA and protein expressions:
for (c in celltype) {
  
  #all_protein_list = get(paste0(c, '_diff_DEGs_DEPs_primary_TN_vs_primary_T'))
  all_protein_list = get(paste0(c, '_diff_DEGs_DEPs_LM_TN_vs_LM_T'))
  
  if(length(all_protein_list) > 0){
    for (p in all_protein_list) {
      
      #plot protein:
      pdf(paste0('Predicted_', c, '_', p, '_LM_TN_vs_LM_T_20230903.pdf'), height = 5, width = 6)
      DefaultAssay(RNA_integrated) <- "ADT"
      p1 = FeaturePlot(RNA_integrated, features = p, min.cutoff = "q05", max.cutoff = "q95", raster = F)
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
      
      if(length(intersect(rownames(RNA_integrated@assays$RNA@counts), gene)) > 0){
        pdf(paste0('RNA_', c, '_', p, '_LM_TN_vs_LM_T_20230903.pdf'), height = 5, width = 6)
        DefaultAssay(RNA_integrated) <- "RNA"
        p1 = FeaturePlot(RNA_integrated, features = gene, min.cutoff = "q05", max.cutoff = "q95", raster = F)
        print(p1)
        dev.off()
      }
      
    }
  }
}



save.image('case_study_metastasis_GSE178318_DEGs_DEPs_20230826.RData')





