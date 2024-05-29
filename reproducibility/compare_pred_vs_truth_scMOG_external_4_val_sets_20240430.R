compare_seurat_pred_vs_true_protein <- function(feature_file = '/Users/ruoqiaochen/Desktop/surface protein prediction/implementation/protein imputation pipeline/test_query_datasets/BM_GSE128639/by_universal_features/gene_cor_mat_tested_GSE128639/predict_denoise_autoencoder/universal_features_gene_cor_mat_GSE128639_unseen_pros_<0.2NA_128dim_encoded.csv',
                                         file_pred_vs_true_format = 'few_shot_protein_specific_add_linear_model_pred_vs_truth_y_unseen_GSE128639_',
                                         file_pred_vs_true_path = '/Users/ruoqiaochen/Desktop/surface protein prediction/implementation/protein imputation pipeline/test_query_datasets/BM_GSE128639/by_universal_features/predict_DNN_unseen_GSE128639/few_shot_learning_protein_specific_model/',
                                         
                                         save_plot_file_path = "/Users/ruoqiaochen/Desktop/surface protein prediction/implementation/protein imputation pipeline/test_query_datasets/BM_GSE128639/by_universal_features/predict_DNN_unseen_GSE128639/few_shot_learning_protein_specific_model/scatter_plot_per_pro/",
                                         save_plot_file_name = '_few_shot_protein_specific_add_linear_model_GSE128639_universal_model_pre_vs_true_universal_features_gene_cor_mat.pdf',
                                         file_cor_name = '/Users/ruoqiaochen/Desktop/surface protein prediction/implementation/protein imputation pipeline/test_query_datasets/BM_GSE128639/by_universal_features/predict_DNN_unseen_GSE128639/few_shot_learning_protein_specific_model/scatter_plot_per_pro/cor_per_pro_universal_few_shot_protein_specific_add_linear_model_GSE128639_universal_features_gene_cor_mat'){
  
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
      xlab('Seurat') +
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
#pancreas_normal_GSM5025059:

compare_seurat_pred_vs_true_protein(file_pred_vs_true_path = '/scMOG/prediction/',
                             
                             file_pred_vs_true_format = '_GSM5025059_y_pred_vs_truth.csv',
                             
                             
                             feature_file = 'seen_protein_names_GSM5025059cell_features_protein_specific_DNN_training_combined_6_sets_20230115.csv',
                             
                             save_plot_file_path = "/scMOG/scatter_plot_per_pro/",
                             
                             save_plot_file_name = '_pancreas_normal_GSM5025059_pre_vs_true_6_combined_training_sets_scMOG_seen_proteins_20240423.pdf',
                             
                             file_cor_name = '/scMOG/scatter_plot_per_pro/cor_per_pro_pancreas_normal_GSM5025059_6_combined_training_sets_scMOG_seen_proteins_20240430')



#pancreas_pancreatitis_GSM5025052:

compare_seurat_pred_vs_true_protein(file_pred_vs_true_path = '/scMOG/prediction/',
                             
                             file_pred_vs_true_format = '_GSM5025052_y_pred_vs_truth.csv',
                             
                             
                             feature_file = 'seen_protein_names_GSM5025052cell_features_protein_specific_DNN_training_combined_6_sets_20230115.csv',
                             
                             save_plot_file_path = "/scMOG/scatter_plot_per_pro/",
                             
                             save_plot_file_name = '_pancreas_pancreatitis_GSM5025052_pre_vs_true_6_combined_training_sets_scMOG_seen_proteins_20240423.pdf',
                             
                             file_cor_name = '/scMOG/scatter_plot_per_pro/cor_per_pro_pancreas_pancreatitis_GSM5025052_6_combined_training_sets_scMOG_seen_proteins_20240430')




#GSE128639:

compare_seurat_pred_vs_true_protein(file_pred_vs_true_path = '/scMOG/prediction/',
                             
                             file_pred_vs_true_format = '_GSE128639_y_pred_vs_truth.csv',
                             

                             feature_file = 'seen_protein_names_GSE128639cell_features_protein_specific_DNN_training_combined_6_sets_20230115.csv',
                             
                             save_plot_file_path = "/scMOG/scatter_plot_per_pro/",
                             
                             save_plot_file_name = '_GSE128639_pre_vs_true_6_combined_training_sets_scMOG_seen_proteins_20240423.pdf',
                             
                             file_cor_name = '/scMOG/scatter_plot_per_pro/cor_per_pro_GSE128639_6_combined_training_sets_scMOG_seen_proteins_20240430')



#lung_COVID_GSE167118:

compare_seurat_pred_vs_true_protein(file_pred_vs_true_path = '/scMOG/prediction/',
                             
                             file_pred_vs_true_format = '_GSE167118_y_pred_vs_truth.csv',
                             
                             
                             feature_file = 'seen_protein_names_GSM5093918cell_features_protein_specific_DNN_training_combined_6_sets_20230115.csv',
                             
                             save_plot_file_path = "/scMOG/scatter_plot_per_pro/",
                             
                             save_plot_file_name = '_BALF_COVID_GSE167118_GSM5093918_pre_vs_true_6_combined_training_sets_scMOG_seen_proteins_20240423.pdf',
                             
                             file_cor_name = '/scMOG/scatter_plot_per_pro/cor_per_pro_BALF_COVID_GSE167118_GSM5093918_6_combined_training_sets_scMOG_seen_proteins_20240430')



