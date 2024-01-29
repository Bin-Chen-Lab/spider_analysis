#! /usr/bin/env Rscript
library(reticulate)
py_config()

library(SAVERX)

#------------------------------------------------------------------------------------------------------------------------
#GSE128639:
data_denoised_QC <- saverx('RNA_QC_count_GSE128639_20230115.csv', ncores = 10) #rownames should be genes, colnames should be cells

data_denoised_estimate_QC <- readRDS(data_denoised_QC)
data_denoised_estimate_QC<- as.data.frame(data_denoised_estimate_QC$estimate)
 

write.csv(data_denoised_estimate_QC, 'RNA_QC_count_GSE128639_denoised_saverx_20230115.csv')

#------------------------------------------------------------------------------------------------------------------------
#GSM5025059, Healthy pancreas:
data_denoised_QC <- saverx('RNA_QC_count_pancreas_normal_GSM5025059_20230115.csv', ncores = 10) #rownames should be genes, colnames should be cells

data_denoised_estimate_QC <- readRDS(data_denoised_QC)
data_denoised_estimate_QC<- as.data.frame(data_denoised_estimate_QC$estimate)


write.csv(data_denoised_estimate_QC, 'RNA_QC_count_GSM5025059_denoised_saverx_20230115.csv')

#------------------------------------------------------------------------------------------------------------------------
#GSM5025052, pancreatic pancreas:
data_denoised_QC <- saverx('RNA_QC_count_pancreas_pancreatitis_GSM5025052_20230115.csv', ncores = 10) #rownames should be genes, colnames should be cells

data_denoised_estimate_QC <- readRDS(data_denoised_QC)
data_denoised_estimate_QC<- as.data.frame(data_denoised_estimate_QC$estimate)


write.csv(data_denoised_estimate_QC, 'pancreas_pancreatitis_GSM5025052/RNA_QC_count_GSM5025052_denoised_saverx_20230115.csv')

#------------------------------------------------------------------------------------------------------------------------
#Annotate lung_COVID_GSE167118:
data_denoised_QC <- saverx('RNA_QC_count_lung_COVID_GSM5093918_20230115.csv', ncores = 10) #rownames should be genes, colnames should be cells

data_denoised_estimate_QC <- readRDS(data_denoised_QC)
data_denoised_estimate_QC<- as.data.frame(data_denoised_estimate_QC$estimate)


write.csv(data_denoised_estimate_QC, 'RNA_QC_count_GSM5093918_denoised_saverx_20230115.csv')

#Rscript SAVER.R &>nohup.out&
#conda activate saverx
#nohup Rscript SAVER.R &

