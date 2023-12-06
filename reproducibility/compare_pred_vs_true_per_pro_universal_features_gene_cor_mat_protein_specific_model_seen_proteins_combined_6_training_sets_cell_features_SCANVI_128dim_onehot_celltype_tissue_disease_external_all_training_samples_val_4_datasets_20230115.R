compare_pred_vs_true_protein <- function(feature_file = '',
                                         file_pred_vs_true_format = '',
                                         file_pred_vs_true_path = '',
                                         
                                         save_plot_file_path = '',
                                         save_plot_file_name = '',
                                         file_cor_name = ''){
  
  library(dplyr)
  library(scales)
  library(ggplot2)
  
  protein_feature <- read.csv(feature_file, stringsAsFactors = FALSE, check.names = F)
  test_set_protein <- protein_feature$protein_name
  
  tmp <- NULL
  spearman_p <- NULL
  RMSE <- NULL
  pearson <- NULL
  pearson_p <- NULL
  r2 <- NULL
  
  for(i in test_set_protein){
    
    file_pred_vs_true <- read.csv(paste0(file_pred_vs_true_path, i, file_pred_vs_true_format), stringsAsFactors = F, row.names = 1)
    
    y_pred <- file_pred_vs_true[, 1]
    y_true <- file_pred_vs_true[, 2]
    
    tmp <- c(tmp, cor.test(y_pred, y_true, method = 'spearman')$estimate)
    spearman_p <-  c(spearman_p, cor.test(y_pred, y_true, method = 'spearman')$p.value)
    RMSE <- c(RMSE, Metrics::rmse(y_true, y_pred))
    pearson <- c(pearson, cor.test(y_pred, y_true, method = 'pearson')$estimate)
    pearson_p <- c(pearson_p, cor.test(y_pred, y_true, method = 'pearson')$p.value)
    r2 <- c(r2, cor(y_pred, y_true)^2)
    
    setwd(save_plot_file_path)
    
    df = data.frame(y_pred = y_pred, y_true = y_true, stringsAsFactors = F)
   
    pdf(paste0(i, '_', length(y_pred), save_plot_file_name), height = 6, width = 6)
    
    plot = ggplot(df, aes(y_pred, y_true)) + 
      geom_point(shape=1, color = 'orange') +
      xlim(0, NA)+
      xlab('SPIDER') +
      ylab('CITE-seq') +
      geom_abline(intercept = 0, slope = 1, color = 'grey', size=1) +
      ggtitle(paste0('cor: ', round(cor.test(y_pred, y_true, method = 'pearson')$estimate, 2))) +
      theme(plot.title = element_text(hjust = 0.5))+ 
      theme_bw()
    
    print(plot)
    dev.off()
    
  }
  
  protein_spearman <- data.frame(spearman = c(tmp, mean(tmp, na.rm = T)), spearman_p = c(spearman_p, mean(spearman_p, na.rm = T)), RMSE = c(RMSE, mean(RMSE, na.rm = T)),  pearson = c(pearson, mean(pearson, na.rm = T)), pearson_p = c(pearson_p, mean(pearson_p, na.rm = T)), r2 = c(r2, mean(r2, na.rm = T)), row.names = c(test_set_protein, 'average'))
  write.csv(protein_spearman, paste0(file_cor_name, '.csv'))
  
  
}

#----------------------------------------------------------------------------------------------------------------------------
#GSE128639:

compare_pred_vs_true_protein(file_pred_vs_true_path = 'predict_DNN_seen_GSE128639/',
                             
                             file_pred_vs_true_format = '_cell_features_protein_specific_onehot_DNN_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_training_combined_6_sets_y_pred_vs_true_seen_GSE128639_universal_features_SCANVI_128dim_impute_HVG_zeros_20230115.csv',
                             
                             
                             feature_file = 'seen_protein_names_GSE128639cell_features_protein_specific_onehot_DNN_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_training_combined_6_sets_SCANVI_128dim_impute_HVG_zeros_20230115.csv',
                             
                             
                             save_plot_file_path = "scatter_plot_per_pro/",
                             
                             save_plot_file_name = '_GSE128639_pre_vs_true_6_combined_training_sets_cell_features_protein_specific_DNN_SCANVI_128dim_impute_HVG_zeros_onehot_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_seen_proteins_20230115.pdf',
                             
                             
                             file_cor_name = 'cor_per_pro_GSE128639_6_combined_training_sets_cell_features_protein_specific_DNN_SCANVI_128dim_impute_HVG_zeros_onehot_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_seen_proteins_20230115')



#pancreas_normal_GSM5025059:

compare_pred_vs_true_protein(file_pred_vs_true_path = 'predict_DNN_seen_pancreas_normal_GSM5025059/',
                             
                             file_pred_vs_true_format = '_cell_features_protein_specific_onehot_DNN_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_training_combined_6_sets_y_pred_vs_true_seen_GSM5025059_universal_features_SCANVI_128dim_impute_HVG_zeros_20230115.csv',
                             
                             
                             feature_file = 'seen_protein_names_GSM5025059cell_features_protein_specific_onehot_DNN_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_training_combined_6_sets_SCANVI_128dim_impute_HVG_zeros_20230115.csv',
                             
                             
                             save_plot_file_path = "scatter_plot_per_pro/",
                             
                             save_plot_file_name = '_pancreas_normal_GSM5025059_pre_vs_true_6_combined_training_sets_cell_features_protein_specific_DNN_SCANVI_128dim_impute_HVG_zeros_onehot_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_seen_proteins_20230115.pdf',
                             
                             
                             file_cor_name = 'cor_per_pro_GSM5025059_6_combined_training_sets_cell_features_protein_specific_DNN_SCANVI_128dim_impute_HVG_zeros_onehot_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_seen_proteins_20230115')



#pancreas_pancreatitis_GSM5025052:

compare_pred_vs_true_protein(file_pred_vs_true_path = 'predict_DNN_seen_GSM5025052_pancreas_pancreatitis/',
                             
                             file_pred_vs_true_format = '_cell_features_protein_specific_onehot_DNN_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_training_combined_6_sets_y_pred_vs_true_seen_GSM5025052_universal_features_SCANVI_128dim_impute_HVG_zeros_20230115.csv',
                             
                             
                             feature_file = 'predict_DNN_seen_GSM5025052_pancreas_pancreatitis/seen_protein_names_GSM5025052cell_features_protein_specific_onehot_DNN_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_training_combined_6_sets_SCANVI_128dim_impute_HVG_zeros_20230115.csv',
                             
                             
                             save_plot_file_path = "scatter_plot_per_pro/",
                             
                             save_plot_file_name = '_pancreas_pancreatitis_GSM5025052_pre_vs_true_6_combined_training_sets_cell_features_protein_specific_DNN_SCANVI_128dim_impute_HVG_zeros_onehot_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_seen_proteins_20230115.pdf',
                             
                             
                             file_cor_name = 'cor_per_pro_GSM5025052_6_combined_training_sets_cell_features_protein_specific_DNN_SCANVI_128dim_impute_HVG_zeros_onehot_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_seen_proteins_20230115')


#lung_COVID_GSE167118:

compare_pred_vs_true_protein(file_pred_vs_true_path = 'predict_DNN_seen_lung_COVID_GSE167118_GSM5093918/',
                             
                             file_pred_vs_true_format = '_cell_features_protein_specific_onehot_DNN_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_training_combined_6_sets_y_pred_vs_true_seen_GSM5093918_universal_features_SCANVI_128dim_impute_HVG_zeros_20230115.csv',
                             
                             
                             feature_file = 'seen_protein_names_GSM5093918cell_features_protein_specific_onehot_DNN_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_training_combined_6_sets_SCANVI_128dim_impute_HVG_zeros_20230115.csv',
                             
                             
                             save_plot_file_path = "scatter_plot_per_pro/",
                             
                             save_plot_file_name = '_lung_COVID_GSE167118_GSM5093918_pre_vs_true_6_combined_training_sets_cell_features_protein_specific_DNN_SCANVI_128dim_impute_HVG_zeros_onehot_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_seen_proteins_20230115.pdf',
                             
                             
                             file_cor_name = 'cor_per_pro_GSM5093918_6_combined_training_sets_cell_features_protein_specific_DNN_SCANVI_128dim_impute_HVG_zeros_onehot_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_seen_proteins_20230115')




