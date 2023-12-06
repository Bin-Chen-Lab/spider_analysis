library(dplyr)
library(Seurat)
#-------------------------------------------------------------
match_training_protein_gene_name = read.csv('all_consistent_720_protein_gene_names_combined_6_training_sets_20230115.csv', stringsAsFactors = F, check.names = F)

match_training_protein_gene_name = match_training_protein_gene_name[, c('consistent_protein_name', 'gene_name')]

match_training_protein_gene_name = unique(match_training_protein_gene_name) #289 proteins
rownames(match_training_protein_gene_name) = match_training_protein_gene_name$consistent_protein_name




featureplot_seurat_RNA_protein <- function(RNA_count_file = 'GSE128639_citeseq_adt_counts.tsv',
                                                   ADT_count_file = 'ADT_QC_CLR_count_GSE128639_20230115.csv',
                                                   dataset_id = 'GSE128639',
                                                   save_plot_file_path = 'RNA_vs_protein/',
                                                   save_plot_file_name = 'RNA_LogNormalized_vs_protein_CLR_truth_20230115.csv',
                                                   model = ...,
                                           file_pred_vs_true_path = 'predict_DNN_seen_pancreas_normal_GSM5025059/',
                                           file_pred_vs_true_format = '_.csv'){
  

  
  if(dataset_id == 'GSE128639'){
    RNA <- as.sparse(read.csv(RNA_count_file, stringsAsFactors = F, row.names = 1, check.names = F,  sep='\t'))
  }
  if(dataset_id == 'GSM5025059'|dataset_id == 'GSM5025052'|dataset_id == 'GSM5093918'){
    RNA <- as.sparse(read.csv(RNA_count_file, stringsAsFactors = F, row.names = 1, check.names = F))
  }
  
  RNA <- CollapseSpeciesExpressionMatrix(RNA)
  RNA <- CreateSeuratObject(counts = RNA)
  
  if(dataset_id == 'GSE128639'|dataset_id == 'GSM5025059'|dataset_id == 'GSM5025052'){
    RNA[["percent.mt"]] <- PercentageFeatureSet(RNA, pattern = '^MT-') #Manually check
  }
  if(dataset_id == 'GSM5093918'){
    RNA[["percent.mt"]] <- PercentageFeatureSet(RNA, pattern = '^MT\\.') #Manually check
  }
  
  RNA <- subset(RNA, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 30)
  features = rownames(RNA)[rowSums(RNA@assays$RNA@counts > 0) >= 5]
  RNA <- subset(RNA, features = features)
  RNA <- NormalizeData(RNA, normalization.method = "LogNormalize", scale.factor = 10000)
  
  ADT_CLR <- read.csv(ADT_count_file, stringsAsFactors = F, row.names = 1, check.names = F)
  ADT_CLR = ADT_CLR[, colnames(RNA)]
  all_protein_list = rownames(ADT_CLR)
  #-------------------------------------------------------------
  #Plot RNA vs protein ground truth:
  if(model == 'RNA'){
  
  ADT_CLR <- as.sparse(ADT_CLR)
  RNA[["ADT"]] <- CreateAssayObject(counts = ADT_CLR)
  RNA@assays$ADT@data = ADT_CLR
  
  RNA <- FindVariableFeatures(RNA)
  RNA <- ScaleData(RNA)
  RNA <- RunPCA(RNA, verbose = FALSE)
  ElbowPlot(RNA, ndims = 40)
  RNA <- FindNeighbors(RNA, dims = 1:40)
  RNA <- FindClusters(RNA, resolution = 0.8)
  RNA <- RunUMAP(RNA, dims = 1:40)
  RNA <- ScaleData(RNA, assay = "ADT")
  
  #Plot all RNA:
  use_all_protein_list <- NULL
  for (n in 1:nrow(ADT_CLR)) {
    if(length(which(rownames(RNA) == match_training_protein_gene_name[rownames(ADT_CLR)[n], 'gene_name'])) > 0){
      use_all_protein_list = c(use_all_protein_list, rownames(ADT_CLR)[n])
    }
  }
  
  pdf(paste0(save_plot_file_path, 'RNA_', dataset_id, '_', save_plot_file_name), height = 6, width = 90)
  DefaultAssay(RNA) <- "RNA"
  p1 = FeaturePlot(RNA, features = match_training_protein_gene_name[use_all_protein_list, 'gene_name'], min.cutoff = "q05", max.cutoff = "q95", ncol = length(use_all_protein_list))
  print(p1)
  dev.off()
  
  pdf(paste0(save_plot_file_path, 'ADT_truth_', dataset_id, '_', save_plot_file_name), height = 6, width = 90)
  DefaultAssay(RNA) <- "ADT"
  p2 = FeaturePlot(RNA, features = use_all_protein_list, min.cutoff = "q05", max.cutoff = "q95", ncol = length(use_all_protein_list))
  #p2 = FeaturePlot(RNA, features = use_all_protein_list, cols = c("lightgrey", "darkgreen"), min.cutoff = "q05", max.cutoff = "q95", ncol = length(use_all_protein_list))
  print(p2)
  dev.off()
  
  }
  #-------------------------------------------------------------
  #Plot SPIDER prediction:
  if(model == 'SPIDER'){
    
    file_pred_vs_true <- NULL
    
    for (i in all_protein_list) {
      
      if(file.exists(paste0(file_pred_vs_true_path, file_pred_vs_true_format, i, '_20230115.csv'))){
    
      tmp = read.csv(paste0(file_pred_vs_true_path, file_pred_vs_true_format, i, '_20230115.csv'), stringsAsFactors = F, row.names = 1)
      
      }
      
      colnames(tmp)[1] = i
      
      if(length(file_pred_vs_true) == 0){
        file_pred_vs_true = tmp
      }else{
        file_pred_vs_true = cbind(file_pred_vs_true, tmp)
      }
      
    }
    
    file_pred_vs_true = t(file_pred_vs_true[, all_protein_list])
    
    RNA[["ADT"]] <- CreateAssayObject(counts = file_pred_vs_true)
    RNA@assays$ADT@data = file_pred_vs_true
    
    RNA <- FindVariableFeatures(RNA)
    RNA <- ScaleData(RNA)
    RNA <- RunPCA(RNA, verbose = FALSE)
    ElbowPlot(RNA, ndims = 40)
    RNA <- FindNeighbors(RNA, dims = 1:40)
    RNA <- FindClusters(RNA, resolution = 0.8)
    RNA <- RunUMAP(RNA, dims = 1:40)
    RNA <- ScaleData(RNA, assay = "ADT")
    #BM_small <- subset(BM_RNA, downsample = 300)
    #BM.adt.markers <- FindAllMarkers(BM_small, assay = "ADT", only.pos = TRUE)
    
    pdf(paste0(save_plot_file_path, 'SPIDER_pred_', dataset_id, '_', save_plot_file_name), height = 6, width = 90)
    DefaultAssay(RNA) <- "ADT"
    p1 = FeaturePlot(RNA, features = all_protein_list, min.cutoff = "q05", max.cutoff = "q95", ncol = length(all_protein_list))
    
    print(p1)
    
    dev.off()
    
  }
  
}

#--------------------------------------------------------------------------------------------------------------
#Run function:

#pancreas_GSM5025052, RNA vs prtein ground truth:
featureplot_seurat_RNA_protein(RNA_count_file = 'pancreas_pancreatitis_BL1_RNA.csv',
                               ADT_count_file = 'ADT_QC_CLR_count_pancreas_pancreatitis_GSM5025052_20230115.csv',
                               dataset_id = 'GSM5025052',
                               save_plot_file_path = 'RNA_vs_protein/',
                               save_plot_file_name = 'featureplot_unseen_proteins_RNA_LogNormalized_vs_protein_CLR_truth_20230115.pdf',
                               model = 'RNA')


#--------------------------------------------------------------------------------------------------------------
#Run function:

#pancreas_GSM5025052, SPIDER predicted:
featureplot_seurat_RNA_protein(RNA_count_file = 'pancreas_pancreatitis_BL1_RNA.csv',
                               ADT_count_file = 'ADT_QC_CLR_count_pancreas_pancreatitis_GSM5025052_20230115.csv',
                               dataset_id = 'GSM5025052',
                               save_plot_file_path = 'figures/',
                               save_plot_file_name = 'featureplot_unseen_proteins_SPIDER_pred_vs_protein_CLR_truth_20230115.pdf',
                               model = 'SPIDER',
                               file_pred_vs_true_path = 'DNN_ensemble/',
                               file_pred_vs_true_format = 'zero_shot_top8_Zi_Zt_cosine_similarity_10189_gene_cell_features_from_GSM5025052_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_SCANVI_128dim_DNN_threshold_0.6_ensemble_pred_vs_truth_y_unseen_')



















