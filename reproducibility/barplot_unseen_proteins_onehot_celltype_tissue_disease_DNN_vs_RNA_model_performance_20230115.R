library(dplyr)
library(ggplot2)
#-------------------------------------------------------------
DNN_file = c("cor_per_pro_zero_shot_top8_Zi_Zt_cosine_similarity_10702_gene_cell_features_from_GSE128639_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_20230115.csv",
             "cor_per_pro_zero_shot_top8_Zi_Zt_cosine_similarity_10373_gene_cell_features_from_GSM5025059_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_20230115.csv",
             "cor_per_pro_zero_shot_top8_Zi_Zt_cosine_similarity_10189_gene_cell_features_from_GSM5025052_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_20230115.csv",
             "cor_per_pro_zero_shot_top8_Zi_Zt_cosine_similarity_10489_gene_cell_features_from_GSM5093918_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_20230115.csv")

RNA_vs_protein_file = c('GSE128639_cor_per_pro_RNA_LogNormalized_vs_protein_CLR_truth_20230115.csv',
                        'GSM5025059_cor_per_pro_RNA_LogNormalized_vs_protein_CLR_truth_20230115.csv',
                        'GSM5025052_cor_per_pro_RNA_LogNormalized_vs_protein_CLR_truth_20230115.csv',
                        'GSM5093918_cor_per_pro_RNA_LogNormalized_vs_protein_CLR_truth_20230115.csv')


#dataset_id = c('GSE128639', 'GSM5025059', 'GSM5025052', 'GSM5093918')
dataset_id = c('healthy bone marrow', 'healthy pancreas', 'pancreatitis pancreas', 'COVID-19 BALF')

#Import proteins' corresponding raw RNA counts:
all_external_datasets_pred_vs_pred_protein_cor = NULL

for (i in 1:length(DNN_file)) {
  
  DNN_pearson_cor <- read.csv(DNN_file[i], stringsAsFactors = F, row.names = 1)
  
  #Set a threshold for maximum dot product, then plot pearson correlation performance:
  threshold_dot_product = 0.85
  #threshold_dot_product = -2
  

  #Check the maximum dot product:
  
  GSE128639_summary_dot_products <- read.csv('max_softmax_Zi_Zt_zero_shot_top_8_10702_gene_cell_features_from_GSE128639_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_unseen_proteins_20230115.csv', stringsAsFactors = F, row.names = 1)
  GSM5025059_summary_dot_products <- read.csv('max_softmax_Zi_Zt_zero_shot_top_8_10373_gene_cell_features_from_GSM5025059_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_unseen_proteins_20230115.csv', stringsAsFactors = F, row.names = 1)
  GSM5025052_summary_dot_products <- read.csv('max_softmax_Zi_Zt_zero_shot_top_8_10189_gene_cell_features_from_GSM5025052_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_unseen_proteins_20230115.csv', stringsAsFactors = F, row.names = 1)
  GSM5093918_summary_dot_products <- read.csv('max_softmax_Zi_Zt_zero_shot_top_8_10489_gene_cell_features_from_GSM5093918_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_unseen_proteins_20230115.csv', stringsAsFactors = F, row.names = 1)
  
  if(i == 1){
  GSE128639_summary_dot_products = GSE128639_summary_dot_products[rownames(DNN_pearson_cor)[1:(nrow(DNN_pearson_cor)-1)], ]
  GSE128639_max_coef_performance <- data.frame(max_coef = as.numeric(GSE128639_summary_dot_products), pred_vs_truth_pearson = DNN_pearson_cor$pearson[1:(nrow(DNN_pearson_cor)-1)], row.names = rownames(DNN_pearson_cor)[1:(nrow(DNN_pearson_cor)-1)], stringsAsFactors = F)
  filter_proteins = GSE128639_max_coef_performance[GSE128639_max_coef_performance$max_coef > threshold_dot_product, ] %>% rownames()
  DNN_pearson_cor = DNN_pearson_cor[filter_proteins, ]
  DNN_pearson_cor['average', ] = apply(DNN_pearson_cor, 2, mean)
  DNN_pearson_cor['median', ] = apply(DNN_pearson_cor[1:(nrow(DNN_pearson_cor) - 1), ], 2, median)
  }
  
  if(i == 2){
  GSM5025059_summary_dot_products = GSM5025059_summary_dot_products[rownames(DNN_pearson_cor)[1:(nrow(DNN_pearson_cor)-1)], ]
  GSM5025059_max_coef_performance <- data.frame(max_coef = as.numeric(GSM5025059_summary_dot_products), pred_vs_truth_pearson = DNN_pearson_cor$pearson[1:(nrow(DNN_pearson_cor)-1)], row.names = rownames(DNN_pearson_cor)[1:(nrow(DNN_pearson_cor)-1)], stringsAsFactors = F)
  filter_proteins = GSM5025059_max_coef_performance[GSM5025059_max_coef_performance$max_coef > threshold_dot_product, ] %>% rownames()
  DNN_pearson_cor = DNN_pearson_cor[filter_proteins, ]
  DNN_pearson_cor['average', ] = apply(DNN_pearson_cor, 2, mean)
  DNN_pearson_cor['median', ] = apply(DNN_pearson_cor[1:(nrow(DNN_pearson_cor) - 1), ], 2, median)
  }
  
  if(i == 3){
  GSM5025052_summary_dot_products = GSM5025052_summary_dot_products[rownames(DNN_pearson_cor)[1:(nrow(DNN_pearson_cor)-1)], ]
  GSM5025052_max_coef_performance <- data.frame(max_coef = as.numeric(GSM5025052_summary_dot_products), pred_vs_truth_pearson = DNN_pearson_cor$pearson[1:(nrow(DNN_pearson_cor)-1)], row.names = rownames(DNN_pearson_cor)[1:(nrow(DNN_pearson_cor)-1)], stringsAsFactors = F)
  filter_proteins = GSM5025052_max_coef_performance[GSM5025052_max_coef_performance$max_coef > threshold_dot_product, ] %>% rownames()
  DNN_pearson_cor = DNN_pearson_cor[filter_proteins, ]
  DNN_pearson_cor['average', ] = apply(DNN_pearson_cor, 2, mean)
  DNN_pearson_cor['median', ] = apply(DNN_pearson_cor[1:(nrow(DNN_pearson_cor) - 1), ], 2, median)
  }
  
  if(i == 4){
  GSM5093918_summary_dot_products = GSM5093918_summary_dot_products[rownames(DNN_pearson_cor)[1:(nrow(DNN_pearson_cor)-1)], ]
  GSM5093918_max_coef_performance <- data.frame(max_coef = as.numeric(GSM5093918_summary_dot_products), pred_vs_truth_pearson = DNN_pearson_cor$pearson[1:(nrow(DNN_pearson_cor)-1)], row.names = rownames(DNN_pearson_cor)[1:(nrow(DNN_pearson_cor)-1)], stringsAsFactors = F)
  filter_proteins = GSM5093918_max_coef_performance[GSM5093918_max_coef_performance$max_coef > threshold_dot_product, ] %>% rownames()
  DNN_pearson_cor = DNN_pearson_cor[filter_proteins, ]
  DNN_pearson_cor['average', ] = apply(DNN_pearson_cor, 2, mean)
  DNN_pearson_cor['median', ] = apply(DNN_pearson_cor[1:(nrow(DNN_pearson_cor) - 1), ], 2, median)
  }
  
  DNN_pearson_cor$group = dataset_id[i]
  DNN_pearson_cor$color = 'SPIDER'
  DNN_pearson_cor$barname = paste0(dataset_id[i], '_SPIDER')
  DNN_pearson_cor = DNN_pearson_cor[, c('pearson', 'RMSE', 'group', 'color', 'barname')]
  
  RNA_vs_protein_pearson_cor <- read.csv(RNA_vs_protein_file[i], stringsAsFactors = F, row.names = 1)
  RNA_vs_protein_pearson_cor$group = dataset_id[i]
  RNA_vs_protein_pearson_cor$color = 'Normalized RNA count'
  RNA_vs_protein_pearson_cor$barname = paste0(dataset_id[i], '_Normalized RNA count')
  RNA_vs_protein_pearson_cor = RNA_vs_protein_pearson_cor[, c('pearson', 'RMSE', 'group', 'color', 'barname')]
  
  
  shared_protein_gene = intersect(rownames(RNA_vs_protein_pearson_cor), rownames(DNN_pearson_cor)) 
  
  if(sum(shared_protein_gene == 'average') > 0){
  shared_protein_gene = shared_protein_gene[-which(shared_protein_gene == 'average')]
  }
  
  
  DNN_pearson_cor = DNN_pearson_cor[shared_protein_gene, ]
  RNA_vs_protein_pearson_cor = RNA_vs_protein_pearson_cor[shared_protein_gene, ]

  truth_vs_pred_protein_cor = rbind(RNA_vs_protein_pearson_cor, DNN_pearson_cor) 
  
  if(length(all_external_datasets_pred_vs_pred_protein_cor) == 0){
    all_external_datasets_pred_vs_pred_protein_cor = truth_vs_pred_protein_cor
  }else{
    all_external_datasets_pred_vs_pred_protein_cor = rbind(all_external_datasets_pred_vs_pred_protein_cor, truth_vs_pred_protein_cor)
  }
  
}


#-------------------------------------------------------------
#Plot pearson correlation:
pdf('RNA_count_vs_model_pred_vs_truth_unseen_protein_cosine_similarities_threshold_0.85_pearson_cor_combined_6_training_sets_all_training_samples_external_4_val_datasets_20230115.pdf', width = 12, height = 8)

library(ggplot2)
library(ggsignif)

plot <- ggplot(all_external_datasets_pred_vs_pred_protein_cor, aes(x=factor(group, levels = c('healthy bone marrow', 'healthy pancreas', 'pancreatitis pancreas', 'COVID-19 BALF')), 
                                                                   y=pearson, 
                                                                   color = factor(color, levels = c('SPIDER', 'Normalized RNA count'))))+
  geom_boxplot(width = 0.5, outlier.shape = NA)+
  theme(legend.position = "right", legend.title = element_blank(), legend.key.size = unit(2, 'cm'))+
  geom_point(aes(x=factor(group, levels = c('healthy bone marrow', 'healthy pancreas', 'pancreatitis pancreas', 'COVID-19 BALF')), 
                 y=pearson, 
                 color = factor(color, levels = c('SPIDER', 'Normalized RNA count'))), 
             size=1, alpha=0.9, position = position_jitterdodge(dodge.width = 0.5)) +
  xlab('')+
  ylab('Correlation')+
  theme(text = element_text(size=20))+ 
  theme_bw()

print(plot)

dev.off()

#-------------------------------------------------------------
#Plot RMSE:
pdf('RNA_count_vs_model_pred_vs_truth_unseen_protein_cosine_similarities_threshold_0.85_RMSE_combined_6_training_sets_all_training_samples_external_4_val_datasets_20230115.pdf', width = 12, height = 8)

library(ggplot2)
library(ggsignif)

plot <- ggplot(all_external_datasets_pred_vs_pred_protein_cor, aes(x=factor(group, levels = c('healthy bone marrow', 'healthy pancreas', 'pancreatitis pancreas', 'COVID-19 BALF')), 
                                                                   y=RMSE, 
                                                                   color = factor(color, levels = c('SPIDER', 'Normalized RNA count'))))+
  geom_boxplot(width = 0.5, outlier.shape = NA)+
  theme(legend.position = "right", legend.title = element_blank(), legend.key.size = unit(2, 'cm'))+
  geom_point(aes(x=factor(group, levels = c('healthy bone marrow', 'healthy pancreas', 'pancreatitis pancreas', 'COVID-19 BALF')), 
                 y=RMSE, 
                 color = factor(color, levels = c('SPIDER', 'Normalized RNA count'))), 
             size=1, alpha=0.9, position = position_jitterdodge(dodge.width = 0.5)) +
  coord_cartesian(y = c(0, 5))+
  xlab('')+
  ylab('RMSE')+
  #coord_cartesian(ylim=c(0, max(filter(all_external_datasets_pred_vs_pred_protein_cor, color != 'Normalized RNA count')$RMSE)+0.2)) +
  theme(text = element_text(size=20))+ 
  theme_bw()

print(plot)

dev.off()








