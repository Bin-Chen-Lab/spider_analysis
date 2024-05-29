import pandas as pd
from sklearn.metrics import mean_squared_error, r2_score
import numpy as np
import time
from sklearn.ensemble import RandomForestRegressor
import pickle
from sklearn.metrics import r2_score
import warnings
import math
import os
import random

#------------------------------------------------------------------------------------------------------------------------------
#prediction:
#changeable parameters:
#GSE128639:
data_denoised_estimate_QC_file = 'RNA_QC_count_GSE128639_denoised_saverx_20230115.csv'
ADT_CLR_file_truth = 'ADT_QC_CLR_count_GSE128639_20230115.csv'
protein_list_file = 'seen_protein_names_GSE128639cell_features_protein_specific_onehot_DNN_celltype_tissue_disease_64_32_16_0.0001_training_combined_6_sets_SCANVI_128dim_20230115.csv'
dataset_id = 'GSE128639'

#------------------------------------------------------------------------------------------------------------------------------
#prediction:
#changeable parameters:
#GSM5025059,  Healthy pancreas:
data_denoised_estimate_QC_file = 'RNA_QC_count_GSM5025059_denoised_saverx_20230115.csv'
ADT_CLR_file_truth = 'ADT_QC_CLR_count_pancreas_normal_GSM5025059_20230115.csv'
dataset_id = 'GSM5025059'
#------------------------------------------------------------------------------------------------------------------------------
#prediction:
#changeable parameters:
#GSM5025052, pancreatitis pancreas:
data_denoised_estimate_QC_file = 'RNA_QC_count_GSM5025052_denoised_saverx_20230115.csv'
ADT_CLR_file_truth = 'ADT_QC_CLR_count_pancreas_pancreatitis_GSM5025052_20230115.csv'
dataset_id = 'GSM5025052'
#------------------------------------------------------------------------------------------------------------------------------
#prediction:
#changeable parameters:
#lung_COVID_GSE167118:
data_denoised_estimate_QC_file = 'RNA_QC_count_GSM5093918_denoised_saverx_20230115.csv'
ADT_CLR_file_truth = 'ADT_QC_CLR_count_lung_COVID_GSM5093918_20230115.csv'
dataset_id = 'GSE167118'
#------------------------------------------------------------------------------------------------------------------------------
RNA_test = pd.read_csv(data_denoised_estimate_QC_file, index_col = 'Unnamed: 0')
#test_protein_list = pd.read_csv(protein_list_file)['protein_name'].tolist()
test_ADT_truth = pd.read_csv(ADT_CLR_file_truth, index_col = 'Unnamed: 0')
test_protein_list = test_ADT_truth.index.tolist()

for p in test_protein_list:
    RNA_protein_corr = pd.read_csv(p + '_RNA_feature_selection_cor_ref_combined_6_datasets_all_training_samples_20230115.csv', index_col = 'Unnamed: 0')
    shared_RNA_features = list(set(RNA_protein_corr.index) & set(RNA_test.index))
    RNA_protein_corr = RNA_protein_corr.loc[shared_RNA_features, :]
    rna_sort = RNA_protein_corr.sort_values('corr', ascending=False)
    high_corr_mrna = rna_sort[0:20].index.tolist()
    x_test = RNA_test.loc[high_corr_mrna, :].T
    y_test =  test_ADT_truth.loc[p, :].T
    with open(p + '_' + dataset_id, 'rb') as f:
        rf = pickle.load(f)
        df = pd.DataFrame({'y_pred': rf.predict(x_test), 'y_truth': y_test})
        df.to_csv(p + '_' + dataset_id + '_ref_combined_6_datasets_all_training_samples_20230115.csv')


