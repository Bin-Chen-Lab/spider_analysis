library(dplyr)
#-------------------------------------------------------------
seurat_file = c("cor_per_pro_GSE128639_6_combined_training_sets_all_training_samples_SeuratV4_seen_proteins_20230115.csv",
                "cor_per_pro_pancreas_normal_GSM5025059_6_combined_training_sets_all_training_samples_SeuratV4_seen_proteins_20230115.csv",
                "cor_per_pro_pancreas_pancreatitis_GSM5025052_6_combined_training_sets_all_training_samples_SeuratV4_seen_proteins_20230115.csv",
                'cor_per_pro_BALF_COVID_GSE167118_GSM5093918_6_combined_training_sets_all_training_samples_SeuratV4_seen_proteins_20230115.csv')


DNN_file = c("cor_per_pro_GSE128639_6_combined_training_sets_cell_features_protein_specific_DNN_SCANVI_128dim_onehot_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_seen_proteins_20230115.csv",
             "cor_per_pro_GSM5025059_6_combined_training_sets_cell_features_protein_specific_DNN_SCANVI_128dim_onehot_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_seen_proteins_20230115.csv",
             "cor_per_pro_GSM5025052_6_combined_training_sets_cell_features_protein_specific_DNN_SCANVI_128dim_onehot_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_seen_proteins_20230115.csv",
             "cor_per_pro_GSM5093918_6_combined_training_sets_cell_features_protein_specific_DNN_SCANVI_128dim_onehot_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_seen_proteins_20230115.csv")

RNA_vs_protein_file = c('GSE128639_cor_per_pro_RNA_LogNormalized_vs_protein_CLR_truth_20230115.csv',
                        'GSM5025059_cor_per_pro_RNA_LogNormalized_vs_protein_CLR_truth_20230115.csv',
                        'GSM5025052_cor_per_pro_RNA_LogNormalized_vs_protein_CLR_truth_20230115.csv',
                        'GSM5093918_cor_per_pro_RNA_LogNormalized_vs_protein_CLR_truth_20230115.csv')

totalvi_file = c("cor_per_pro_GSE128639_combined_6_training_sets_all_training_samples_ep400_TOTALVI_seen_proteins_20230115.csv",
                "cor_per_pro_GSM5025059_combined_6_training_sets_all_training_samples_ep400_TOTALVI_seen_proteins_20230115.csv",
                "cor_per_pro_GSM5025052_combined_6_training_sets_all_training_samples_ep400_TOTALVI_seen_proteins_20230115.csv",
                'cor_per_pro_GSM5093918_combined_6_training_sets_all_training_samples_ep400_TOTALVI_seen_proteins_20230115.csv')

imputed_RNA_vs_protein_file = c('GSE128639_cor_per_pro_saverx_imputed_RNA_vs_protein_CLR_truth_20230912.csv',
                        'GSM5025059_cor_per_pro_saverx_imputed_RNA_vs_protein_CLR_truth_20230912.csv',
                        'GSM5025052_cor_per_pro_saverx_imputed_RNA_vs_protein_CLR_truth_20230912.csv',
                        'GSM5093918_cor_per_pro_saverx_imputed_RNA_vs_protein_CLR_truth_20230912.csv')

sciPENN_file = c("cor_per_pro_GSE128639_combined_6_training_sets_sciPENN_seen_proteins_20231126.csv",
                 "cor_per_pro_GSM5025059_combined_6_training_sets_all_training_samples_sciPENN_seen_proteins_20231126.csv",
                 "cor_per_pro_GSM5025052_combined_6_training_sets_all_training_samples_sciPENN_seen_proteins_20231126.csv",
                 'cor_per_pro_GSM5093918_combined_6_training_sets_sciPENN_seen_proteins_20231126.csv')

#dataset_id = c('GSE128639', 'GSM5025059', 'GSM5025052', 'GSM5093918')
dataset_id = c('healthy bone marrow', 'healthy pancreas', 'pancreatitis pancreas', 'COVID-19 BALF')

#Import proteins' corresponding raw RNA counts:
all_external_datasets_pred_vs_pred_protein_cor = NULL

for (i in 1:length(seurat_file)) {
  
  seurat_pearson_cor <- read.csv(seurat_file[i], stringsAsFactors = F, row.names = 1)
  seurat_pearson_cor$group = dataset_id[i]
  seurat_pearson_cor$color = 'Seurat'
  seurat_pearson_cor$barname = paste0(dataset_id[i], '_Seurat')
  seurat_pearson_cor = seurat_pearson_cor[, c('pearson', 'RMSE', 'group', 'color', 'barname')]
  
  DNN_pearson_cor <- read.csv(DNN_file[i], stringsAsFactors = F, row.names = 1)
  DNN_pearson_cor$group = dataset_id[i]
  DNN_pearson_cor$color = 'SPIDER'
  DNN_pearson_cor$barname = paste0(dataset_id[i], '_SPIDER')
  DNN_pearson_cor = DNN_pearson_cor[, c('pearson', 'RMSE', 'group', 'color', 'barname')]
  
  RNA_vs_protein_pearson_cor <- read.csv(RNA_vs_protein_file[i], stringsAsFactors = F, row.names = 1)
  RNA_vs_protein_pearson_cor$group = dataset_id[i]
  RNA_vs_protein_pearson_cor$color = 'Normalized RNA count'
  RNA_vs_protein_pearson_cor$barname = paste0(dataset_id[i], '_Normalized RNA count')
  RNA_vs_protein_pearson_cor = RNA_vs_protein_pearson_cor[, c('pearson', 'RMSE', 'group', 'color', 'barname')]
  
  totalvi_pearson_cor <- read.csv(totalvi_file[i], stringsAsFactors = F, row.names = 1)
  totalvi_pearson_cor$group = dataset_id[i]
  totalvi_pearson_cor$color = 'totalVI'
  totalvi_pearson_cor$barname = paste0(dataset_id[i], '_totalVI')
  totalvi_pearson_cor = totalvi_pearson_cor[, c('pearson', 'RMSE', 'group', 'color', 'barname')]
  
  imputed_RNA_vs_protein_pearson_cor <- read.csv(imputed_RNA_vs_protein_file[i], stringsAsFactors = F, row.names = 1)
  imputed_RNA_vs_protein_pearson_cor$group = dataset_id[i]
  imputed_RNA_vs_protein_pearson_cor$color = 'Imputed RNA'
  imputed_RNA_vs_protein_pearson_cor$barname = paste0(dataset_id[i], '_Imputed RNA')
  imputed_RNA_vs_protein_pearson_cor = imputed_RNA_vs_protein_pearson_cor[, c('pearson', 'RMSE', 'group', 'color', 'barname')]
  
  sciPENN_pearson_cor <- read.csv(sciPENN_file[i], stringsAsFactors = F, row.names = 1)
  sciPENN_pearson_cor$group = dataset_id[i]
  sciPENN_pearson_cor$color = 'sciPENN'
  sciPENN_pearson_cor$barname = paste0(dataset_id[i], '_sciPENN')
  sciPENN_pearson_cor = sciPENN_pearson_cor[, c('pearson', 'RMSE', 'group', 'color', 'barname')]
  
  
  shared_protein_gene = intersect(rownames(seurat_pearson_cor), rownames(DNN_pearson_cor)) %>% intersect(rownames(RNA_vs_protein_pearson_cor)) %>% intersect(rownames(imputed_RNA_vs_protein_pearson_cor)) %>% intersect(rownames(sciPENN_pearson_cor))
  
  if(sum(shared_protein_gene == 'average') > 0){
  shared_protein_gene = shared_protein_gene[-which(shared_protein_gene == 'average')]
  }
  
  seurat_pearson_cor = seurat_pearson_cor[shared_protein_gene, ]
  DNN_pearson_cor = DNN_pearson_cor[shared_protein_gene, ]
  RNA_vs_protein_pearson_cor = RNA_vs_protein_pearson_cor[shared_protein_gene, ]
  totalvi_pearson_cor = totalvi_pearson_cor[shared_protein_gene, ]
  imputed_RNA_vs_protein_pearson_cor = imputed_RNA_vs_protein_pearson_cor[shared_protein_gene, ]
  sciPENN_pearson_cor = sciPENN_pearson_cor[shared_protein_gene, ]
  
  truth_vs_pred_protein_cor = rbind(seurat_pearson_cor, DNN_pearson_cor) %>% rbind(RNA_vs_protein_pearson_cor) %>% rbind(totalvi_pearson_cor) %>% rbind(imputed_RNA_vs_protein_pearson_cor) %>% rbind(sciPENN_pearson_cor)
  
  if(length(all_external_datasets_pred_vs_pred_protein_cor) == 0){
    all_external_datasets_pred_vs_pred_protein_cor = truth_vs_pred_protein_cor
  }else{
    all_external_datasets_pred_vs_pred_protein_cor = rbind(all_external_datasets_pred_vs_pred_protein_cor, truth_vs_pred_protein_cor)
  }
  
}

#-------------------------------------------------------------
#Plot pearson correlation:
pdf('RNA_count_vs_model_vs_baseline_pred_vs_truth_seen_protein_pearson_cor_combined_6_training_sets_all_training_samples_external_4_val_datasets_20230115.pdf', width = 16, height = 8)

library(ggplot2)
library(ggsignif)

plot <- ggplot(all_external_datasets_pred_vs_pred_protein_cor, aes(x=factor(group, levels = c('healthy bone marrow', 'healthy pancreas', 'pancreatitis pancreas', 'COVID-19 BALF')), 
                                                                   y=pearson, 
                                                                   color = factor(color, levels = c('SPIDER', 'Seurat', 'totalVI', 'sciPENN', 'Normalized RNA count', 'Imputed RNA'))))+
  geom_boxplot(width = 0.5, outlier.shape = NA)+
  theme(legend.position = "right", legend.title = element_blank(), legend.key.size = unit(2, 'cm'))+
  geom_point(aes(x=factor(group, levels = c('healthy bone marrow', 'healthy pancreas', 'pancreatitis pancreas', 'COVID-19 BALF')), 
                 y=pearson, 
                 color = factor(color, levels = c('SPIDER', 'Seurat', 'totalVI', 'sciPENN', 'Normalized RNA count', 'Imputed RNA'))), 
                 size=1, alpha=0.9, position = position_jitterdodge(dodge.width = 0.5)) +
  xlab('')+
  ylab('Correlation')+
  theme(text = element_text(size=20))+ 
  theme_bw()

print(plot)

dev.off()

#-------------------------------------------------------------
#Plot RMSE:
pdf('RNA_count_vs_model_vs_baseline_pred_vs_truth_seen_protein_RMSE_combined_6_training_sets_all_training_samples_external_4_val_datasets_20230115.pdf', width = 16, height = 8)

library(ggplot2)
library(ggsignif)

plot <- ggplot(all_external_datasets_pred_vs_pred_protein_cor, aes(x=factor(group, levels = c('healthy bone marrow', 'healthy pancreas', 'pancreatitis pancreas', 'COVID-19 BALF')), 
                                                                   y=RMSE, 
                                                                   color = factor(color, levels = c('SPIDER', 'Seurat', 'totalVI', 'sciPENN', 'Normalized RNA count', 'Imputed RNA'))))+
  geom_boxplot(width = 0.5, outlier.shape = NA)+
  theme(legend.position = "right", legend.title = element_blank(), legend.key.size = unit(2, 'cm'))+
  geom_point(aes(x=factor(group, levels = c('healthy bone marrow', 'healthy pancreas', 'pancreatitis pancreas', 'COVID-19 BALF')), 
                 y=RMSE, 
                 color = factor(color, levels = c('SPIDER', 'Seurat', 'totalVI', 'sciPENN', 'Normalized RNA count', 'Imputed RNA'))), 
                 size=1, alpha=0.9, position = position_jitterdodge(dodge.width = 0.5)) +
  coord_cartesian(y = c(0, 5))+
  xlab('')+
  ylab('RMSE')+
  theme(text = element_text(size=20))+ 
  theme_bw()

print(plot)

dev.off()






