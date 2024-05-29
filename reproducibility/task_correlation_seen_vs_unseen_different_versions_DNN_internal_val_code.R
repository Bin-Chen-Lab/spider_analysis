library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)
#----------------------------------------------------------------------------------------------------------------------------
DNN_V1_internal_performance_pearson <- read.csv('cor_per_pro_internal_val_6_combined_training_sets_cell_features_protein_specific_DNN_SCANVI_10dim_128_32_64_0.0001_seen_proteins_20230115.csv', stringsAsFactors = F, row.names = 1)
DNN_V2_internal_performance_pearson <- read.csv('cor_per_pro_internal_val_6_combined_training_sets_cell_features_protein_specific_DNN_onehot_celltype_tissue_disease_SCANVI_128dim_64_32_16_0.0001_seen_proteins_20230115.csv', stringsAsFactors = F, row.names = 1)
DNN_V3_internal_performance_pearson <- read.csv('cor_per_pro_internal_val_6_combined_training_sets_cell_features_protein_specific_DNN_SCANVI_128dim_64_32_16_0.0001_seen_proteins_20230115.csv', stringsAsFactors = F, row.names = 1)
DNN_V4_internal_performance_pearson <- read.csv('cor_per_pro_internal_val_6_combined_training_sets_cell_features_protein_specific_DNN_SCANVI_256dim_64_32_16_0.0001_seen_proteins_20230115.csv', stringsAsFactors = F, row.names = 1)
DNN_V5_internal_performance_pearson <- read.csv('cor_per_pro_internal_val_6_combined_training_sets_cell_features_protein_specific_DNN_SCANVI_64dim_64_32_16_0.0001_seen_proteins_20230115.csv', stringsAsFactors = F, row.names = 1)
DNN_V6_internal_performance_pearson <- read.csv('cor_per_pro_internal_val_6_combined_training_sets_cell_features_protein_specific_DNN_onehot_celltype_SCANVI_128dim_64_32_16_0.0001_seen_proteins_20230115.csv', stringsAsFactors = F, row.names = 1)
DNN_V7_internal_performance_pearson <- read.csv('cor_per_pro_internal_val_6_combined_training_sets_cell_features_protein_specific_DNN_onehot_tissue_SCANVI_128dim_64_32_16_0.0001_seen_proteins_20230115.csv', stringsAsFactors = F, row.names = 1)
DNN_V8_internal_performance_pearson <- read.csv('cor_per_pro_internal_val_6_combined_training_sets_cell_features_protein_specific_DNN_onehot_disease_SCANVI_128dim_64_32_16_0.0001_seen_proteins_20230115.csv', stringsAsFactors = F, row.names = 1)

DNN_V1_internal_performance_pearson$model_id = 'Model1' 
DNN_V5_internal_performance_pearson$model_id = 'Model2' 
DNN_V4_internal_performance_pearson$model_id = 'Model3'
DNN_V3_internal_performance_pearson$model_id = 'Model4' 
DNN_V2_internal_performance_pearson$model_id = 'Model5' 
DNN_V6_internal_performance_pearson$model_id = 'Model6' 
DNN_V7_internal_performance_pearson$model_id = 'Model7' 
DNN_V8_internal_performance_pearson$model_id = 'Model8' 

DNN_internal_performance_pearson_seen = c(median(DNN_V1_internal_performance_pearson[-nrow(DNN_V1_internal_performance_pearson), ]$pearson, na.rm = T), 
                                         median(DNN_V5_internal_performance_pearson[-nrow(DNN_V5_internal_performance_pearson), ]$pearson, na.rm = T), 
                                         median(DNN_V4_internal_performance_pearson[-nrow(DNN_V4_internal_performance_pearson), ]$pearson, na.rm = T), 
                                     median(DNN_V3_internal_performance_pearson[-nrow(DNN_V3_internal_performance_pearson), ]$pearson, na.rm = T), 
                                     median(DNN_V2_internal_performance_pearson[-nrow(DNN_V2_internal_performance_pearson), ]$pearson, na.rm = T), 
                                     median(DNN_V6_internal_performance_pearson[-nrow(DNN_V6_internal_performance_pearson), ]$pearson, na.rm = T), 
                                     median(DNN_V7_internal_performance_pearson[-nrow(DNN_V7_internal_performance_pearson), ]$pearson, na.rm = T),
                                     median(DNN_V8_internal_performance_pearson[-nrow(DNN_V8_internal_performance_pearson), ]$pearson, na.rm = T)
                                         )

DNN_V1_internal_performance_pearson_unseen <- read.csv('cor_per_pro_universal_zero_shot_top8_Zi_Zt_cosine_similarity_10702_gene_training_combined_6_sets_64_32_16_0.0001_SCANVI_10dim_DNN_internal_val_threshold_0.6_ensemble_internal_all_samples_pre_vs_true_unseen_20240430.csv', stringsAsFactors = F, row.names = 1)
DNN_V2_internal_performance_pearson_unseen <- read.csv('cor_per_pro_universal_zero_shot_top8_Zi_Zt_cosine_similarity_10702_gene_training_combined_6_sets_64_32_16_0.0001_SCANVI_64dim_DNN_internal_val_threshold_0.6_ensemble_internal_all_samples_pre_vs_true_unseen_20240430.csv', stringsAsFactors = F, row.names = 1)
DNN_V3_internal_performance_pearson_unseen <- read.csv('cor_per_pro_universal_zero_shot_top8_Zi_Zt_cosine_similarity_10702_gene_training_combined_6_sets_64_32_16_0.0001_SCANVI_256dim_DNN_internal_val_threshold_0.6_ensemble_internal_all_samples_pre_vs_true_unseen_20240430.csv', stringsAsFactors = F, row.names = 1)
DNN_V4_internal_performance_pearson_unseen <- read.csv('cor_per_pro_universal_zero_shot_top8_Zi_Zt_cosine_similarity_10702_gene_training_combined_6_sets_64_32_16_0.0001_SCANVI_128dim_DNN_internal_val_threshold_0.6_ensemble_internal_all_samples_pre_vs_true_unseen_20240430.csv', stringsAsFactors = F, row.names = 1)
DNN_V5_internal_performance_pearson_unseen <- read.csv('cor_per_pro_universal_zero_shot_top8_Zi_Zt_cosine_similarity_10702_gene_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_internal_val_threshold_0.6_ensemble_internal_all_samples_pre_vs_true_unseen_20230115.csv', stringsAsFactors = F, row.names = 1)
DNN_V6_internal_performance_pearson_unseen <- read.csv('cor_per_pro_universal_zero_shot_top8_Zi_Zt_cosine_similarity_10702_gene_training_combined_6_sets_64_32_16_0.0001_SCANVI_128dim_onehot_celltype_DNN_internal_val_threshold_0.6_ensemble_internal_all_samples_pre_vs_true_unseen_20240430.csv', stringsAsFactors = F, row.names = 1)
DNN_V7_internal_performance_pearson_unseen <- read.csv('cor_per_pro_universal_zero_shot_top8_Zi_Zt_cosine_similarity_10702_gene_training_combined_6_sets_64_32_16_0.0001_SCANVI_128dim_onehot_tissue_DNN_internal_val_threshold_0.6_ensemble_internal_all_samples_pre_vs_true_unseen_20240430.csv', stringsAsFactors = F, row.names = 1)
DNN_V8_internal_performance_pearson_unseen <- read.csv('cor_per_pro_universal_zero_shot_top8_Zi_Zt_cosine_similarity_10702_gene_training_combined_6_sets_64_32_16_0.0001_SCANVI_128dim_onehot_disease_DNN_internal_val_threshold_0.6_ensemble_internal_all_samples_pre_vs_true_unseen_20240430.csv', stringsAsFactors = F, row.names = 1)

DNN_internal_performance_pearson_unseen = c(mean(DNN_V1_internal_performance_pearson_unseen[-nrow(DNN_V1_internal_performance_pearson_unseen), ]$pearson, na.rm = T), 
                                          mean(DNN_V2_internal_performance_pearson_unseen[-nrow(DNN_V2_internal_performance_pearson_unseen), ]$pearson, na.rm = T), 
                                          mean(DNN_V3_internal_performance_pearson_unseen[-nrow(DNN_V3_internal_performance_pearson_unseen), ]$pearson, na.rm = T), 
                                          mean(DNN_V4_internal_performance_pearson_unseen[-nrow(DNN_V4_internal_performance_pearson_unseen), ]$pearson, na.rm = T), 
                                          mean(DNN_V5_internal_performance_pearson_unseen[-nrow(DNN_V5_internal_performance_pearson_unseen), ]$pearson, na.rm = T), 
                                          mean(DNN_V6_internal_performance_pearson_unseen[-nrow(DNN_V6_internal_performance_pearson_unseen), ]$pearson, na.rm = T), 
                                          mean(DNN_V7_internal_performance_pearson_unseen[-nrow(DNN_V7_internal_performance_pearson_unseen), ]$pearson, na.rm = T),
                                          mean(DNN_V8_internal_performance_pearson_unseen[-nrow(DNN_V8_internal_performance_pearson_unseen), ]$pearson, na.rm = T)
)
#----------------------------------------------------------------------------------------------------------------------------
#Plot task correlations:
pdf('task_correlation_seen_vs_unseen_internal_val_performance_20240430.pdf', width = 8, height = 8)

library(ggplot2)
library(ggsignif)

p <- plot(DNN_internal_performance_pearson_seen, DNN_internal_performance_pearson_unseen)

print(p)

dev.off()


