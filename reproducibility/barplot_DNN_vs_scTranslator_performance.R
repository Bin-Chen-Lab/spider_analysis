#-------------------------------------------------------------
#change parameter:
DNN_file = c("cor_per_pro_GSE128639_6_combined_training_sets_cell_features_protein_specific_DNN_SCANVI_128dim_onehot_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_seen_proteins_20230115.csv",
             "cor_per_pro_GSM5025059_6_combined_training_sets_cell_features_protein_specific_DNN_SCANVI_128dim_onehot_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_seen_proteins_20230115.csv",
             "cor_per_pro_GSM5025052_6_combined_training_sets_cell_features_protein_specific_DNN_SCANVI_128dim_onehot_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_seen_proteins_20230115.csv",
             "cor_per_pro_GSM5093918_6_combined_training_sets_cell_features_protein_specific_DNN_SCANVI_128dim_onehot_celltype_tissue_disease_all_training_samples_64_32_16_0.0001_seen_proteins_20230115.csv")


sctranslator_file = c("cor_per_pro_GSE128639_pretrained_scTranslator_proteins_20230115.csv",
                "cor_per_pro_GSM5025059_pretrained_scTranslator_proteins_20230115.csv",
                "cor_per_pro_GSM5025052_pretrained_scTranslator_proteins_20230115.csv",
                'cor_per_pro_GSM5093918_pretrained_scTranslator_proteins_20230115.csv')


#-------------------------------------------------------------
#dataset_id = c('GSE128639', 'GSM5025059', 'GSM5025052', 'GSM5093918')
dataset_id = c('GSE128639', 'GSM5025059', 'GSM5025052', 'GSM5093918')

#Import proteins' corresponding raw RNA counts:
all_external_datasets_pred_vs_pred_protein_cor = NULL

for (i in 1:length(DNN_file)) {
  
  DNN_pearson_cor <- read.csv(DNN_file[i], stringsAsFactors = F, row.names = 1)
  DNN_pearson_cor$group = dataset_id[i]
  DNN_pearson_cor$color = 'SPIDER'
  DNN_pearson_cor$barname = paste0(dataset_id[i], '_SPIDER')
  DNN_pearson_cor = DNN_pearson_cor[, c('pearson', 'RMSE', 'group', 'color', 'barname')]
  
  sctranslator_pearson_cor <- read.csv(sctranslator_file[i], stringsAsFactors = F, row.names = 1)
  sctranslator_pearson_cor$group = dataset_id[i]
  sctranslator_pearson_cor$color = 'scTranslator'
  sctranslator_pearson_cor$barname = paste0(dataset_id[i], '_scTranslator')
  sctranslator_pearson_cor = sctranslator_pearson_cor[, c('pearson', 'RMSE', 'group', 'color', 'barname')]
  
  
  shared_protein_gene = intersect(rownames(sctranslator_pearson_cor), rownames(DNN_pearson_cor)) 
  
  if(sum(shared_protein_gene == 'average') > 0){
  shared_protein_gene = shared_protein_gene[-which(shared_protein_gene == 'average')]
  }
  
  DNN_pearson_cor = DNN_pearson_cor[shared_protein_gene, ]
  sctranslator_pearson_cor = sctranslator_pearson_cor[shared_protein_gene, ]
  
  truth_vs_pred_protein_cor = rbind(sctranslator_pearson_cor, DNN_pearson_cor) 
  
  if(length(all_external_datasets_pred_vs_pred_protein_cor) == 0){
    all_external_datasets_pred_vs_pred_protein_cor = truth_vs_pred_protein_cor
  }else{
    all_external_datasets_pred_vs_pred_protein_cor = rbind(all_external_datasets_pred_vs_pred_protein_cor, truth_vs_pred_protein_cor)
  }
  
}


#-------------------------------------------------------------
#Plot pearson correlation:
pdf('scTranslator_vs_model_pred_vs_truth_seen_protein_pearson_cor_external_4_val_datasets_20230506.pdf', width = 12, height = 8)

library(ggplot2)
library(ggsignif)

plot <- ggplot(all_external_datasets_pred_vs_pred_protein_cor, aes(x=factor(group), y=pearson, color = factor(color)))+
  geom_boxplot(width = 0.5, outlier.shape = NA)+
  theme(legend.position = "right", legend.title = element_blank(), legend.key.size = unit(2, 'cm'))+
  #geom_jitter(size=1, alpha=0.9, position = position_jitterdodge(dodge.width = 0)) +
  geom_point(aes(x=factor(group), y=pearson, color = factor(color)), size=2, alpha=0.9, position = position_jitterdodge(dodge.width = 0.5)) +
  xlab('External validation set')+
  ylab('Pearson correlation')+
  theme(text = element_text(size=20))+
  theme_bw()

print(plot)

dev.off()

#-------------------------------------------------------------
#Plot RMSE:
pdf('scTranslator_vs_model_pred_vs_truth_seen_protein_RMSE_external_4_val_datasets_20230506.pdf', width = 12, height = 8)

library(ggplot2)
library(ggsignif)

plot <- ggplot(all_external_datasets_pred_vs_pred_protein_cor, aes(x=factor(group), y=RMSE, color = factor(color)))+
  geom_boxplot(width = 0.5, outlier.shape = NA)+
  theme(legend.position = "right", legend.title = element_blank(), legend.key.size = unit(2, 'cm'))+
  geom_point(aes(x=factor(group), y=RMSE, color = factor(color)), size=2, alpha=0.9, position = position_jitterdodge(dodge.width = 0.5)) +
  xlab('External validation set')+
  ylab('RMSE')+
  #coord_cartesian(ylim=c(0, max(filter(all_external_datasets_pred_vs_pred_protein_cor, color != 'Normalized RNA count')$RMSE)+0.2)) +
  theme(text = element_text(size=20))+ 
  theme_bw()

print(plot)

dev.off()






