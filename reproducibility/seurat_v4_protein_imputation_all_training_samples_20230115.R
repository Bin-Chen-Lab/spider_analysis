#Seurat v4 Reference Mapping:
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)

#Import RNA & ADT data for combined all 6 training sets:
load('CLR_ADT_&_RNA_combined_6_training_datasets_20230115.RData')

ref_RNA = NormalizeData(RNA, normalization.method = "LogNormalize", scale.factor = 10000)
ref_RNA <- FindVariableFeatures(ref_RNA)

#-----------------------------------------------------------------------------------------------------
#changeable parameters:
#Import query data set: GSM5025059:
RNA_count_file = 'pancreas_normal_BL8_RNA.csv'
RNA <- as.sparse(read.csv(RNA_count_file, stringsAsFactors = F, row.names = 1, check.names = F))
RNA <- CollapseSpeciesExpressionMatrix(RNA)
RNA <- CreateSeuratObject(counts = RNA)
RNA[["percent.mt"]] <- PercentageFeatureSet(RNA, pattern = '^MT-') #Manually check
RNA <- subset(RNA, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 30)
features = rownames(RNA)[rowSums(RNA@assays$RNA@counts > 0) >= 5]
RNA <- subset(RNA, features = features)
RNA <- NormalizeData(RNA, normalization.method = "LogNormalize", scale.factor = 10000)
RNA <- FindVariableFeatures(RNA)
RNA <- ScaleData(RNA)
RNA <- RunPCA(RNA, verbose = FALSE)
RNA <- FindNeighbors(RNA, dims = 1:40)
RNA <- FindClusters(RNA, resolution = 0.8)
RNA <- RunUMAP(RNA, dims = 1:40)
all_protein_list <- c('CD3','CD19','CD25','CD56','CD206','CD11c','CD45','CD11b','CD197','CD45RA','HLA-DR','CD4','CD8')
test_dataset_id = 'GSM5025059'
save_pred_vs_truth_path = '/compare_seen_protein_baseline_model/Seurat v4/'

DNN_file_pred_vs_true_path = '/predict_DNN_seen_pancreas_normal_GSM5025059/'
DNN_file_pred_vs_true_format = '_cell_features_protein_specific_DNN_training_combined_6_sets_y_pred_vs_true_seen_GSM5025059_universal_features_gene_cor_mat_20230115.csv'

#-----------------------------------------------------------------------------------------------------
#changeable parameters:
#Import query data set: GSM5025052:
RNA_count_file = 'pancreas_pancreatitis_BL1_RNA.csv'
RNA <- as.sparse(read.csv(RNA_count_file, stringsAsFactors = F, row.names = 1, check.names = F))
RNA <- CollapseSpeciesExpressionMatrix(RNA)
RNA <- CreateSeuratObject(counts = RNA)
RNA[["percent.mt"]] <- PercentageFeatureSet(RNA, pattern = '^MT-') #Manually check
RNA <- subset(RNA, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 30)
features = rownames(RNA)[rowSums(RNA@assays$RNA@counts > 0) >= 5]
RNA <- subset(RNA, features = features)
RNA <- NormalizeData(RNA, normalization.method = "LogNormalize", scale.factor = 10000)
RNA <- FindVariableFeatures(RNA)
RNA <- ScaleData(RNA)
RNA <- RunPCA(RNA, verbose = FALSE)
RNA <- FindNeighbors(RNA, dims = 1:40)
RNA <- FindClusters(RNA, resolution = 0.8)
RNA <- RunUMAP(RNA, dims = 1:40)
all_protein_list <- c('CD3','CD19','CD25','CD56','CD206','CD11c','CD45','CD11b','CD197','CD45RA','HLA-DR','CD4','CD8')
test_dataset_id = 'GSM5025052'
save_pred_vs_truth_path = '/compare_seen_protein_baseline_model/Seurat v4/'

DNN_file_pred_vs_true_path = '/predict_DNN_seen_GSM5025052_pancreas_pancreatitis/'
DNN_file_pred_vs_true_format = '_cell_features_protein_specific_onehot_DNN_celltype_tissue_disease_64_32_16_0.0001_training_combined_6_sets_y_pred_vs_true_seen_GSM5025052_universal_features_SCANVI_128dim_20230115.csv'
#-----------------------------------------------------------------------------------------------------
#changeable parameters:
#Import query data set: GSE128639:
RNA_count_file = 'GSE128639_citeseq_rna_counts.tsv'
RNA <- as.sparse(read.csv(RNA_count_file, stringsAsFactors = F, row.names = 1, check.names = F,  sep='\t'))
RNA <- CollapseSpeciesExpressionMatrix(RNA)
RNA <- CreateSeuratObject(counts = RNA)
RNA[["percent.mt"]] <- PercentageFeatureSet(RNA, pattern = '^MT-') #Manually check
RNA <- subset(RNA, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 30)
features = rownames(RNA)[rowSums(RNA@assays$RNA@counts > 0) >= 5]
RNA <- subset(RNA, features = features)
RNA <- NormalizeData(RNA, normalization.method = "LogNormalize", scale.factor = 10000)
RNA <- FindVariableFeatures(RNA)
RNA <- ScaleData(RNA)
RNA <- RunPCA(RNA, verbose = FALSE)
RNA <- FindNeighbors(RNA, dims = 1:40)
RNA <- FindClusters(RNA, resolution = 0.8)
RNA <- RunUMAP(RNA, dims = 1:40)
all_protein_list <- read.csv('seen_protein_names_GSE128639cell_features_protein_specific_onehot_DNN_celltype_tissue_disease_64_32_16_0.0001_training_combined_6_sets_SCANVI_128dim_20230115.csv', stringsAsFactors = FALSE, check.names = F)$protein_name
test_dataset_id = 'GSE128639'
save_pred_vs_truth_path = '/compare_seen_protein_baseline_model/Seurat v4/'

DNN_file_pred_vs_true_path = '/predict_DNN_seen_GSE128639/'
DNN_file_pred_vs_true_format = '_cell_features_protein_specific_DNN_training_combined_6_sets_y_pred_vs_true_seen_GSE128639_universal_features_gene_cor_mat_20230115.csv'
#-----------------------------------------------------------------------------------------------------
#changeable parameters:
#Import query data set: GSM5093918:
RNA_count_file = 'lung_RNA_Sar01L_CD3pos.csv'
RNA <- as.sparse(read.csv(RNA_count_file, stringsAsFactors = F, row.names = 1, check.names = F))
RNA <- CollapseSpeciesExpressionMatrix(RNA)
RNA <- CreateSeuratObject(counts = RNA)
RNA[["percent.mt"]] <- PercentageFeatureSet(RNA, pattern = '^MT\\.') #Manually check
RNA <- subset(RNA, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 30)
features = rownames(RNA)[rowSums(RNA@assays$RNA@counts > 0) >= 5]
RNA <- subset(RNA, features = features)
RNA <- NormalizeData(RNA, normalization.method = "LogNormalize", scale.factor = 10000)
RNA <- FindVariableFeatures(RNA)
RNA <- ScaleData(RNA)
RNA <- RunPCA(RNA, verbose = FALSE)
RNA <- FindNeighbors(RNA, dims = 1:40)
RNA <- FindClusters(RNA, resolution = 0.8)
RNA <- RunUMAP(RNA, dims = 1:40)
all_protein_list <- read.csv('seen_protein_names_GSM5093918cell_features_protein_specific_onehot_DNN_celltype_tissue_disease_64_32_16_0.0001_training_combined_6_sets_SCANVI_128dim_20230115.csv', stringsAsFactors = FALSE, check.names = F)$protein_name
test_dataset_id = 'GSM5093918'
save_pred_vs_truth_path = '/compare_seen_protein_baseline_model/Seurat v4/'

DNN_file_pred_vs_true_path = '/predict_DNN_seen_lung_COVID_GSE167118_GSM5093918/'
DNN_file_pred_vs_true_format = '_cell_features_protein_specific_DNN_training_combined_6_sets_y_pred_vs_true_seen_GSM5093918_universal_features_gene_cor_mat_20230115.csv'


#-----------------------------------------------------------------------------------------------------

#Import reference data set:
ref_RNA <- ScaleData(ref_RNA)
ref_RNA <- RunPCA(ref_RNA, verbose = FALSE)
ref_RNA <- FindNeighbors(ref_RNA, dims = 1:40)
ref_RNA <- FindClusters(ref_RNA, resolution = 0.8)
ref_RNA <- RunUMAP(ref_RNA, dims = 1:40)

for(p in all_protein_list){ #One protein one model
  
  train_cell_id <- read.csv(paste0(p, '_cell_features_training_combined_6_sets_all_training_samples_cell_id_20230115.csv'), stringsAsFactors = F)$train_cell_id #Import training set cell id
  reference = ref_RNA[, train_cell_id]
  
  ADT_normalized <- NULL
  
  for(j in c('files_13165_BM_all_ADT', 'files_BM_AML_GSE143363_all_ADT', 'files_peritoneum_GSE172155_all_ADT',
             'files_pleura_GSE172155_all_ADT', 'files_GBM_all_ADT', 'files_PBMC_all_ADT')){
    
    if(sum(rownames(get(j)) == p) > 0){
      if(length(ADT_normalized) != 0){
        ADT_normalized = cbind(ADT_normalized, as.data.frame(get(j)[p, ]))
      }else{
        ADT_normalized = as.data.frame(get(j)[p, ])
      }
    }
    
    }#Import ADT data
  
  ADT_normalized = ADT_normalized[, train_cell_id]
    
  reference[["ADT"]] <- CreateAssayObject(counts = ADT_normalized) 
  reference@assays$ADT@data <- as.matrix(ADT_normalized) 


DefaultAssay(RNA) <- "RNA"
DefaultAssay(reference) <- "RNA"


anchors <- FindTransferAnchors(
  reference = reference,
  query = RNA,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:40
)

RNA_query <- MapQuery(
  anchorset = anchors,
  query = RNA,
  reference = reference,
  refdata = list(
    predicted_ADT = "ADT"
  ),
  reference.reduction = "pca")
  #reduction.model = "umap")


ADT_pred <- t(as.data.frame(RNA_query@assays$predicted_ADT@data)) #1 protein

DNN_file_pred_vs_true <- read.csv(paste0(DNN_file_pred_vs_true_path, p, DNN_file_pred_vs_true_format), stringsAsFactors = F, row.names = 1)
seurat_pred_vs_true = data.frame(y_pred = ADT_pred[, p], y_truth = DNN_file_pred_vs_true[rownames(ADT_pred), ]$y_truth, stringsAsFactors = F)

write.csv(seurat_pred_vs_true, paste0(save_pred_vs_truth_path, 'SeuratV4_pred_vs_truth_', p, '_', test_dataset_id, '_ref_combined_6_datasets_all_training_samples_20230115.csv'))

}



