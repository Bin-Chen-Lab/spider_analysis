library(Seurat)
library(dplyr)
library(ggplot2)
library(enrichR)

load('case_study_metastasis_GSE178318_DEGs_DEPs_20230826.RData')

#Find DEGs / DEPs between different sites (tissues) using all cells:
group1 = 'primary_TN'
group2 = 'LM_TN'

#
#Use only highly confident proteins:
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

seen_all_protein_list <- read.csv('protein_gene_names_union_289_DNNs_from_combined_6_training_sets_DNN_onehot_celltype_tissue_disease_SCANVI_128dim_internal_val_threshold_0.6_20230115.csv', stringsAsFactors = F)
use_unseen_proteins = setdiff(rownames(unseen_file_pred), seen_all_protein_list$gene_name) 
use_unseen_proteins = setdiff(use_unseen_proteins, seen_all_protein_list$consistent_protein_name) #excluding overlap from 289 seen proteins
use_unseen_proteins = intersect(use_unseen_proteins[get(paste0('unseen_protein_max_similarity_', group1))[use_unseen_proteins, 'max_inferred_coef'] > 0.85], use_unseen_proteins[get(paste0('unseen_protein_max_similarity_', group2))[use_unseen_proteins, 'max_inferred_coef'] > 0.85])

#
#For enrichment pathway analysis, all the protein names should be changed to gene names:
#Check duplicated gene names for multiple proteins:
seen_file_pred <- seen_file_pred[use_seen_proteins, ]
unseen_file_pred <- unseen_file_pred[use_unseen_proteins, ]

match_training_protein_gene_name = read.csv('protein_gene_names_union_289_DNNs_from_combined_6_training_sets_DNN_onehot_celltype_tissue_disease_SCANVI_128dim_internal_val_threshold_0.6_20230115.csv', stringsAsFactors = F)
rownames(match_training_protein_gene_name) = match_training_protein_gene_name$consistent_protein_name


new_row_names <- rownames(seen_file_pred)

for (n in 1:nrow(seen_file_pred)) {
  
  if(length(intersect(rownames(match_training_protein_gene_name), rownames(seen_file_pred)[n])) > 0){
    if(match_training_protein_gene_name[rownames(seen_file_pred)[n] , 'gene_name'] != ' '){
      new_row_names[n] = match_training_protein_gene_name[rownames(seen_file_pred)[n] , 'gene_name']
    }
  }
  
}

rownames(seen_file_pred)[which(duplicated(new_row_names))] #Manually check
seen_file_pred = seen_file_pred[setdiff(rownames(seen_file_pred), c("CD45", "CD45RO", "CD45RB", "HLA-ABC", "Folate")), ] #Manually check
#
for (n in 1:nrow(seen_file_pred)) {
  
  if(length(intersect(rownames(match_training_protein_gene_name), rownames(seen_file_pred)[n])) > 0){
    if(match_training_protein_gene_name[rownames(seen_file_pred)[n] , 'gene_name'] != ' '){
      rownames(seen_file_pred)[n] = match_training_protein_gene_name[rownames(seen_file_pred)[n] , 'gene_name']
    }
  }
  
}

which(duplicated(rownames(seen_file_pred), rownames(unseen_file_pred)))
seen_unseeen_file_pred <- rbind(seen_file_pred, unseen_file_pred) #304 highly confident proteins

shared_protein_RNA <- intersect(rownames(seen_unseeen_file_pred), rownames(RNA_integrated[['RNA']])) #282 highly confident proteins with corresponding genes

DefaultAssay(RNA_integrated) = "RNA"
RNA2 <- subset(RNA_integrated, features = shared_protein_RNA)

RNA2[["ADT"]] <- CreateAssayObject(counts = seen_unseeen_file_pred[shared_protein_RNA, ])
RNA2@assays$ADT@data = seen_unseeen_file_pred[shared_protein_RNA, ]
RNA2 <- ScaleData(RNA2, assay = "ADT")

#
#Plot enrichment pathways:
group1 = 'primary_treatment_naïve'
group2 = 'liver_metastases_treatment_naïve'
#
#Based on predicted ADT:
celltype_group1 = as.data.frame(table(RNA2@meta.data$cell_type[RNA2@meta.data$study == 'primary_treatment_naïve']))
celltype_group1 <- as.character(celltype_group1$Var1[celltype_group1$Freq >=3])
celltype_group2 <- as.data.frame(table(RNA2@meta.data$cell_type[RNA2@meta.data$study == 'liver_metastases_treatment_naïve']))
celltype_group2 <- as.character(celltype_group2$Var1[celltype_group2$Freq >=3])

celltype = intersect(celltype_group1, celltype_group2)

for(c in celltype){
  
pdf(paste0('enrichment_pathway_ADT_BP_primary_TN_vs_LM_TN_', c, '_20230905.pdf'), height = 8, width = 16)
p1 <- DEenrichRPlot(
  RNA2,
  ident.1 = colnames(RNA2)[RNA2@meta.data$study == group2 & RNA2@meta.data$cell_type == c],
  ident.2 = colnames(RNA2)[RNA2@meta.data$study == group1 & RNA2@meta.data$cell_type == c],
  balanced = FALSE,
  assay = 'ADT',
  enrich.database = c('GO_Biological_Process_2023'),
  #enrich.database = c('GO_Cellular_Component_2023'),
  #enrich.database = c('GO_Molecular_Function_2023'),
  max.genes = 1000,
  return.gene.list = F
  #cols = 'indianred2'
)
print(p1)
dev.off()

Sys.sleep(2)

}

#
#Based on RNA, using highly-confident proteins' corresponding genes:
for(c in celltype){
  
pdf(paste0('enrichment_pathway_RNA_BP_primary_TN_vs_LM_TN_', c, '_20230905.pdf'), height = 8, width = 16)
p2 <- DEenrichRPlot(
  RNA2,
  ident.1 = colnames(RNA2)[RNA2@meta.data$study == group2 & RNA2@meta.data$cell_type == c],
  ident.2 = colnames(RNA2)[RNA2@meta.data$study == group1 & RNA2@meta.data$cell_type == c],
  balanced = FALSE,
  assay = 'RNA',
  enrich.database = c('GO_Biological_Process_2023'),
  #enrich.database = c('GO_Cellular_Component_2023'),
  #enrich.database = c('GO_Molecular_Function_2023'),
  max.genes = 1000,
  return.gene.list = F
)
print(p2)
dev.off()

Sys.sleep(2)

}

#
#Based on RNA, using all genes:
for(c in celltype){
  
  pdf(paste0('enrichment_pathway_use_all_genes_RNA_BP_primary_TN_vs_LM_TN_', c, '_20230905.pdf'), height = 8, width = 20)
  p2 <- DEenrichRPlot(
    RNA_integrated,
    ident.1 = colnames(RNA_integrated)[RNA_integrated@meta.data$study == group1 & RNA_integrated@meta.data$cell_type == c],
    ident.2 = colnames(RNA_integrated)[RNA_integrated@meta.data$study == group2 & RNA_integrated@meta.data$cell_type == c],
    balanced = FALSE,
    assay = 'RNA',
    enrich.database = c('GO_Biological_Process_2023'),
    #enrich.database = c('GO_Cellular_Component_2023'),
    #enrich.database = c('GO_Molecular_Function_2023'),
    max.genes = 1000,
    return.gene.list = F
  )
  print(p2)
  dev.off()
  
  Sys.sleep(2)
  
}

#--------------------------------------------------------------------------------------------------------------------
#Alluvial plot:
use_celltype <- c("CAFs", "Epithelial cells", "pDCs")

all_pos.er <- NULL

#positive markers:
for (c in use_celltype) {
  
object = RNA2
#object = RNA_integrated
ident.1 = colnames(RNA2)[RNA2@meta.data$study == group2 & RNA2@meta.data$cell_type == c]
ident.2 = colnames(RNA2)[RNA2@meta.data$study == group1 & RNA2@meta.data$cell_type == c]
balanced = T
#assay = 'RNA'
assay = 'ADT'
enrich.database = c('GO_Biological_Process_2023')
max.genes = 1000
return.gene.list = F
logfc.threshold = 0.25
test.use = 'wilcox'
p.val.cutoff = 0.05
cols = NULL
num.pathway = 5


  enrichr.installed <- PackageCheck("enrichR", error = FALSE)

  assay <- assay %||% DefaultAssay(object = object)
  
  DefaultAssay(object = object) <- assay
  
  all.markers <- FindMarkers(
    object = object,
    ident.1 = ident.1,
    ident.2 = ident.2,
    only.pos = FALSE,
    logfc.threshold = logfc.threshold,
    test.use = test.use,
    assay = assay
  )
  
  pos.markers <- all.markers[all.markers[, 2] > logfc.threshold & all.markers[, 1] < p.val.cutoff, , drop = FALSE]
  
  if(nrow(pos.markers) == 0){
    message("No positive markers pass the logfc.thershold")
    pos.er <- c()
  }else{
    pos.markers.list <- rownames(x = pos.markers)[1:min(max.genes, nrow(x = pos.markers))]
    pos.er <- enrichR::enrichr(genes = pos.markers.list, databases = enrich.database)
    pos.er <- do.call(what = cbind, args = pos.er)
    pos.er$log10pval <- -log10(x = pos.er[, paste(enrich.database, sep = ".", "P.value")])
    pos.er$term <- pos.er[, paste(enrich.database, sep = ".", "Term")]
    pos.er <- pos.er[1:num.pathway, ]
    pos.er$term <- factor(x = pos.er$term, levels = pos.er$term[order(pos.er$log10pval)])
    gene.list <- list(pos = pos.er)
    
    pos.er$cell_type = c
    
    if(length(all_pos.er) == 0){
      all_pos.er = pos.er
    }else{
      all_pos.er = rbind(all_pos.er, pos.er)
    }
    
  }
  Sys.sleep(2)
}

#negative markers:
all_neg.er <- NULL

for (c in use_celltype) {
  
  
  object = RNA2
  #object = RNA_integrated
  ident.1 = colnames(RNA2)[RNA2@meta.data$study == group2 & RNA2@meta.data$cell_type == c]
  ident.2 = colnames(RNA2)[RNA2@meta.data$study == group1 & RNA2@meta.data$cell_type == c]
   balanced = T
  #assay = 'RNA'
  assay = 'ADT'
  enrich.database = c('GO_Biological_Process_2023')
  max.genes = 1000
  return.gene.list = F
  logfc.threshold = 0.25
  test.use = 'wilcox'
  p.val.cutoff = 0.05
  cols = NULL
  num.pathway = 5
  
  
  enrichr.installed <- PackageCheck("enrichR", error = FALSE)
  
  assay <- assay %||% DefaultAssay(object = object)
  
  DefaultAssay(object = object) <- assay
  
  all.markers <- FindMarkers(
    object = object,
    ident.1 = ident.1,
    ident.2 = ident.2,
    only.pos = FALSE,
    logfc.threshold = logfc.threshold,
    test.use = test.use,
    assay = assay
  )
  
  neg.markers <- all.markers[all.markers[, 2] < logfc.threshold & all.markers[, 1] < p.val.cutoff, , drop = FALSE]
  neg.markers.list <- rownames(x = neg.markers)[1:min(max.genes, nrow(x = neg.markers))]
  neg.er <- enrichR::enrichr(genes = neg.markers.list, databases = enrich.database)
  neg.er <- do.call(what = cbind, args = neg.er)
  neg.er$log10pval <- -log10(x = neg.er[, paste(enrich.database, sep = ".", "P.value")])
  neg.er$term <- neg.er[, paste(enrich.database, sep = ".", "Term")]
  neg.er <- neg.er[1:num.pathway, ]
  neg.er$term <- factor(x = neg.er$term, levels = neg.er$term[order(neg.er$log10pval)])
  
  neg.er$cell_type = c
  if(length(all_neg.er) == 0){
    all_neg.er = neg.er
  }else{
    all_neg.er = rbind(all_neg.er, neg.er)
  }
  
  Sys.sleep(2)
}


#Alluvial plot:
library(ggalluvial)

pdf('enrichment_pathway_predicted_BP_primary_TN_vs_LM_TN_all_celltypes_positive_DEPs_20230913.pdf', height = 80, width = 100)
p1 = ggplot(data = all_pos.er,
       aes(axis1 = factor(cell_type), axis2 = factor(term), y = log10pval)) +
  geom_alluvium(aes(fill = factor(term))) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)), size=8) +
  scale_x_discrete(limits = c("cell_type", "term")) +
  theme_void()+
  theme(legend.position = "none")
print(p1)
dev.off()


pdf('enrichment_pathway_predicted_BP_primary_TN_vs_LM_TN_all_celltypes_negative_DEPs_20230913.pdf', height = 80, width = 100)
p2 = ggplot(data = all_neg.er,
            aes(axis1 = factor(cell_type), axis2 = factor(term), y = log10pval)) +
  geom_alluvium(aes(fill = factor(term))) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)), size=8) +
  scale_x_discrete(limits = c("cell_type", "term")) +
  theme_void()+
  theme(legend.position = "none")
print(p2)
dev.off()






  
  
  
  
 
