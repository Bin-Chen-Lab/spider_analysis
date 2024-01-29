#! /usr/bin/env Rscript
library(reticulate)
py_config()

library(cTPnet)
library(Seurat)

model_file_path="/compare_seen_protein_baseline_model/cTPnet/cTPnet_default_pretrained/cTPnet_weight_24"
data_type='dataframe'

#------------------------------------------------------------------------------------------------------------------------
#GSE128639:
data_denoised_estimate_QC <- read.csv('RNA_QC_count_GSE128639_denoised_saverx_20230115.csv', stringsAsFactors = F, row.names = 1)
#y_pred = CreateSeuratObject(data_denoised_estimate_QC)
y_pred = cTPnet(data_denoised_estimate_QC, data_type, model_file_path)
rownames(y_pred) = unlist(strsplit(rownames(y_pred), 'ctpnet_'))[seq(2, 48, 2)]
rownames(y_pred)[15] = 'CD127'
rownames(y_pred)[18] = 'CD278'
rownames(y_pred)[24] = 'HLA-DR'

save(y_pred, file = 'cTPnet_default_pretrained_24_proteins_pred_y_GSE128639.RData')

#------------------------------------------------------------------------------------------------------------------------
#GSM5025059, Healthy pancreas:
data_denoised_estimate_QC <- read.csv('RNA_QC_count_GSM5025059_denoised_saverx_20230115.csv', stringsAsFactors = F, row.names = 1)
y_pred = cTPnet(data_denoised_estimate_QC, data_type, model_file_path)
rownames(y_pred) = unlist(strsplit(rownames(y_pred), 'ctpnet_'))[seq(2, 48, 2)]
rownames(y_pred)[15] = 'CD127'
rownames(y_pred)[18] = 'CD278'
rownames(y_pred)[24] = 'HLA-DR'

save(y_pred, file = 'cTPnet_default_pretrained_24_proteins_pred_y_GSM5025059.RData')

#------------------------------------------------------------------------------------------------------------------------
#GSM5025052, pancreatic pancreas:
data_denoised_estimate_QC <- read.csv('pancreas_pancreatitis_GSM5025052/RNA_QC_count_GSM5025052_denoised_saverx_20230115.csv', stringsAsFactors = F, row.names = 1)
y_pred = cTPnet(data_denoised_estimate_QC, data_type, model_file_path)
rownames(y_pred) = unlist(strsplit(rownames(y_pred), 'ctpnet_'))[seq(2, 48, 2)]
rownames(y_pred)[15] = 'CD127'
rownames(y_pred)[18] = 'CD278'
rownames(y_pred)[24] = 'HLA-DR'

save(y_pred, file = 'cTPnet_default_pretrained_24_proteins_pred_y_GSM5025052.RData')

#------------------------------------------------------------------------------------------------------------------------
#Annotate lung_COVID_GSE167118:
data_denoised_estimate_QC <- read.csv('RNA_QC_count_GSM5093918_denoised_saverx_20230115.csv', stringsAsFactors = F, row.names = 1)
y_pred = cTPnet(data_denoised_estimate_QC, data_type, model_file_path)
rownames(y_pred) = unlist(strsplit(rownames(y_pred), 'ctpnet_'))[seq(2, 48, 2)]
rownames(y_pred)[15] = 'CD127'
rownames(y_pred)[18] = 'CD278'
rownames(y_pred)[24] = 'HLA-DR'

save(y_pred, file = 'cTPnet_default_pretrained_24_proteins_pred_y_GSM5093918.RData')


#Rscript SAVER.R &>nohup.out&
#conda activate saverx
#nohup Rscript SAVER.R &

