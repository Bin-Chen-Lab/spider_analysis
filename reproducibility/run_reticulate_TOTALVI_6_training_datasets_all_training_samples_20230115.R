library(Seurat)
library(dplyr)
library(anndata)

library(reticulate)
sc <- import('scanpy', convert = FALSE)
scvi <- import('scvi', convert = FALSE)
setwd('/compare_seen_protein_baseline_model/TOTALVI/')
plt <- import('matplotlib.pyplot', convert = FALSE)
np <- import('numpy', convert = FALSE)
os <- import('os', convert = FALSE)
anndata <- import('anndata', convert = FALSE)
pd <- import('pandas', convert = FALSE)
time <- import('time', convert = FALSE)
torch <- import('torch', convert = FALSE)


scvi$settings$progress_bar_style = 'tqdm'

#Import RNA & ADT data for combined all 6 training sets:
load('CLR_ADT_&_RNA_combined_6_training_datasets_20230115.RData')
ref_RNA = RNA

#Combine all PBMC raw ADT files
files_ADT = NULL

for (i in c('P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8')) {
  for (j in c(0, 2, 7)) {
    files_ADT = c(files_ADT, paste0('PBMC_ADT_', i, '_', j))
    }
}

files_RNA = ls()[grep('GSM5008738_seuratobj', ls())]

files_PBMC_all_raw_ADT = get(files_ADT[1])[, colnames(get(files_RNA[1]))]
for (i in 2:length(files_ADT)) {
  files_PBMC_all_raw_ADT = cbind(files_PBMC_all_raw_ADT, get(files_ADT[i])[, colnames(get(files_RNA[i]))])
}

colnames(files_PBMC_all_raw_ADT) = colnames(files_PBMC_all_ADT)
rownames(files_PBMC_all_raw_ADT) = rownames(files_PBMC_all_ADT)

colnames(BM_ADT) = colnames(files_13165_BM_all_ADT)
rownames(BM_ADT) = rownames(files_13165_BM_all_ADT)

colnames(GBM_ADT) = colnames(files_GBM_all_ADT)
rownames(GBM_ADT) = rownames(files_GBM_all_ADT)

colnames(peritoneum_ADT) = colnames(files_peritoneum_GSE172155_all_ADT)
rownames(peritoneum_ADT) = rownames(files_peritoneum_GSE172155_all_ADT)

colnames(pleura_ADT) = colnames(files_pleura_GSE172155_all_ADT)
rownames(pleura_ADT) = rownames(files_pleura_GSE172155_all_ADT)

colnames(BM_AML_ADT) = colnames(files_BM_AML_GSE143363_all_ADT)
rownames(BM_AML_ADT) = rownames(files_BM_AML_GSE143363_all_ADT)

#-------------------------------------------------------------------------------------------------
#changeable parameters:
#Import query data set GSM5025059:
save_file_path = '/compare_seen_protein_baseline_model/TOTALVI/'
test_dataset_id = 'GSM5025059'
#all_protein_list <- c('CD3','CD19','CD25','CD56','CD206','CD11c','CD45','CD11b','CD197','CD45RA','HLA-DR','CD4','CD8')
all_protein_list <- c('CD4','CD8')


RNA_count_file = 'pancreas_normal_BL8_RNA.csv'
RNA <- as.sparse(read.csv(RNA_count_file, stringsAsFactors = F, row.names = 1, check.names = F))
RNA <- CollapseSpeciesExpressionMatrix(RNA)
RNA <- CreateSeuratObject(counts = RNA)
RNA[["percent.mt"]] <- PercentageFeatureSet(RNA, pattern = '^MT-') #Manually check
RNA <- subset(RNA, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 30)
features = rownames(RNA)[rowSums(RNA@assays$RNA@counts > 0) >= 5]
RNA <- subset(RNA, features = features)
RNA[['study']] = test_dataset_id

#-------------------------------------------------------------------------------------------------
#changeable parameters:
#Import query data set: GSM5025052:
save_file_path = '/compare_seen_protein_baseline_model/TOTALVI/'
test_dataset_id = 'GSM5025052'
#all_protein_list <- c('CD3','CD19','CD25','CD56','CD206','CD11c','CD45','CD11b','CD197','CD45RA','HLA-DR','CD4','CD8')
all_protein_list <- c('CD19', 'CD56','CD206','CD11c','CD45','CD11b','CD197','CD45RA','HLA-DR','CD4','CD8')

RNA_count_file = 'pancreas_pancreatitis_BL1_RNA.csv'
RNA <- as.sparse(read.csv(RNA_count_file, stringsAsFactors = F, row.names = 1, check.names = F))
RNA <- CollapseSpeciesExpressionMatrix(RNA)
RNA <- CreateSeuratObject(counts = RNA)
RNA[["percent.mt"]] <- PercentageFeatureSet(RNA, pattern = '^MT-') #Manually check
RNA <- subset(RNA, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 30)
features = rownames(RNA)[rowSums(RNA@assays$RNA@counts > 0) >= 5]
RNA <- subset(RNA, features = features)
RNA[['study']] = test_dataset_id

#-------------------------------------------------------------------------------------------------
#changeable parameters:
#Import query data set: GSE128639:
save_file_path = '/compare_seen_protein_baseline_model/TOTALVI/'
test_dataset_id = 'GSE128639'
all_protein_list <- read.csv('seen_protein_names_GSE128639cell_features_protein_specific_onehot_DNN_celltype_tissue_disease_64_32_16_0.0001_training_combined_6_sets_SCANVI_128dim_20230115.csv', stringsAsFactors = FALSE, check.names = F)$protein_name[23:25]


RNA_count_file = 'GSE128639_citeseq_rna_counts.tsv'
RNA <- as.sparse(read.csv(RNA_count_file, stringsAsFactors = F, row.names = 1, check.names = F,  sep='\t'))
RNA <- CollapseSpeciesExpressionMatrix(RNA)
RNA <- CreateSeuratObject(counts = RNA)
RNA[["percent.mt"]] <- PercentageFeatureSet(RNA, pattern = '^MT-') #Manually check
RNA <- subset(RNA, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 30)
features = rownames(RNA)[rowSums(RNA@assays$RNA@counts > 0) >= 5]
RNA <- subset(RNA, features = features)
RNA[['study']] = test_dataset_id

#-------------------------------------------------------------------------------------------------
#changeable parameters:
#Import query data set: GSM5093918:
save_file_path = '/compare_seen_protein_baseline_model/TOTALVI/'
test_dataset_id = 'GSM5093918'
all_protein_list <- read.csv('seen_protein_names_GSM5093918cell_features_protein_specific_onehot_DNN_celltype_tissue_disease_64_32_16_0.0001_training_combined_6_sets_SCANVI_128dim_20230115.csv', stringsAsFactors = FALSE, check.names = F)$protein_name[29:39]


RNA_count_file = 'lung_RNA_Sar01L_CD3pos.csv'
RNA <- as.sparse(read.csv(RNA_count_file, stringsAsFactors = F, row.names = 1, check.names = F))
RNA <- CollapseSpeciesExpressionMatrix(RNA)
RNA <- CreateSeuratObject(counts = RNA)
RNA[["percent.mt"]] <- PercentageFeatureSet(RNA, pattern = '^MT\\.') #Manually check
RNA <- subset(RNA, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 30)
features = rownames(RNA)[rowSums(RNA@assays$RNA@counts > 0) >= 5]
RNA <- subset(RNA, features = features)
RNA[['study']] = test_dataset_id


#-------------------------------------------------------------------------------------------------
#Import reference data set:
for(p in all_protein_list){ #One protein one model
  
  train_cell_id <- read.csv(paste0(p, '_cell_features_training_combined_6_sets_all_training_samples_cell_id_20230115.csv'), stringsAsFactors = F)$train_cell_id #Import training set cell id
  reference = ref_RNA[, train_cell_id]
  
  ADT_raw <- NULL
  
  for(j in c('files_PBMC_all_raw_ADT', 'BM_ADT', 'GBM_ADT', 'peritoneum_ADT', 'pleura_ADT', 'BM_AML_ADT')){
    
    if(sum(rownames(get(j)) == p) > 0){
      if(length(ADT_raw) != 0){
        ADT_raw = cbind(ADT_raw, as.data.frame(get(j)[p, ]))
      }else{
        ADT_raw = as.data.frame(get(j)[p, ])
      }
    }
    
  }#Import ADT data
  
  ADT_raw = ADT_raw[, train_cell_id]
  
  #reference[["ADT"]] <- CreateAssayObject(counts = ADT_raw) 
  
  adata_ref <- sc$AnnData(
    X   = t(as.matrix(GetAssayData(reference,slot='counts',assay = "RNA"))), #scVI requires raw counts
    obs = reference[[]],
    obsm = list(protein_expression = pd$DataFrame(columns=list(rownames(ADT_raw)), index=colnames(reference), data = t(ADT_raw))),
    var = GetAssay(reference)[[]]
  )
  print(adata_ref) #CD3: 113759 × 10702
  
  
  #Concatenate reference with query data set:
  pro_exp = adata_ref$obsm$`__getitem__`('protein_expression')
  data = np$zeros((c(ncol(RNA), pro_exp$shape[1]))) #CD3: 1239 x 1
  
  adata_query <- sc$AnnData(
    X   = t(as.matrix(GetAssayData(RNA,slot='counts',assay = "RNA"))), #scVI requires raw counts
    obs = RNA[[]],
    obsm = list(protein_expression = pd$DataFrame(columns=list(rownames(ADT_raw)), index=colnames(RNA), data = data)),
    var = GetAssay(RNA)[[]]
  )
  print(adata_query) #CD3: 1239 × 12868
  
  adata_full = anndata$concat(list(adata_ref, adata_query)) #114998 × 10373
  print(adata_full)
  
  ref_cell_id = NULL
  
  adata_ref = adata_full[1:dim(reference)[2]-1]$copy()
  adata_query = adata_full[dim(reference)[2]:(dim(reference)[2]+dim(RNA)[2]-1)]$copy()
  
  sc$pp$highly_variable_genes(
    adata_full,
    n_top_genes=as.integer(4000),
    flavor="seurat_v3",
    batch_key="study",
    subset=TRUE
  )
  
  scvi$data$setup_anndata(adata_full, batch_key="study", protein_expression_obsm_key="protein_expression")
  
  model = scvi$model$TOTALVI(
    adata_full,
    latent_distribution="normal"
  )
  
  torch$cuda$synchronize() # wait for warm-up to finish
  start_epoch = time$time()

  
  max_epochs=400
  #model$train(max_epochs=as.integer(max_epochs))
  model$train()
  
  
  torch$cuda$synchronize()
  end_epoch = time$time()
  time_train = py_to_r(end_epoch) - py_to_r(start_epoch)
  write.csv(time_train, paste0(save_file_path, 'training_time_', p, '_ep', max_epochs, '_TOTALVI_ref_combined_6_datasets_external_val_', test_dataset_id, '20230115.csv'))
  

  #Save training loss:
  plt$plot(model$history["elbo_train"], label="train")
  plt$plot(model$history["elbo_validation"], label="validation")
  plt$title("Negative ELBO over training epochs")
  #plt$ylim(1200, 1400)
  plt$legend()
  plt$savefig(paste0(save_file_path, p, '_train_loss_ep', max_epochs, '_TOTALVI_ref_combined_6_datasets_external_val_', test_dataset_id, '20230115.pdf'))
  plt$close()
  
  model$save(paste0(save_file_path, p, '_trained_model_weights_ep', max_epochs, '_TOTALVI_ref_combined_6_datasets_external_val_', test_dataset_id, '_20230115'))

  train_batch_name = as.character(as.data.frame(table(py_to_r(adata_full$obs)$study))$Var1)
  train_batch_name = train_batch_name[-which(train_batch_name == test_dataset_id)]
  train_cell_id = py_to_r(adata_full$obs)
  test_cell_id = rownames(train_cell_id[train_cell_id$study == test_dataset_id, ])
  
  protein_pred = py_to_r(model$get_normalized_expression(transform_batch = train_batch_name, return_mean=TRUE)[1])
  protein_pred = as.data.frame(protein_pred[test_cell_id, ])
  colnames(protein_pred) = 'y_pred'
  rownames(protein_pred) = colnames(RNA)
  
  X_totalVI = py_to_r(model$get_latent_representation())
  rownames(X_totalVI) = py_to_r(adata_full$obs$index)
  write.csv(X_totalVI, paste0(save_file_path, p, '_X_latent_representations_ep', max_epochs, '_TOTALVI_ref_combined_6_datasets_external_val_', test_dataset_id, '20230115.csv'))
  
  write.csv(protein_pred, paste0(save_file_path, 'pred_vs_truth_', p, '_ep', max_epochs, '_TOTALVI_ref_combined_6_datasets_external_val_', test_dataset_id, '20230115.csv'))
  
  
}







