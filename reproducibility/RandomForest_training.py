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

#Construct the input, note that the selected top 20 RNA features should exist in the corresponding external validation set:
#------------------------------------------------------------------------------------------------------------------------------
#prediction:
#changeable parameters:
#GSE128639:
data_denoised_estimate_QC_file = 'RNA_QC_count_GSE128639_denoised_saverx_20230115.csv'
ADT_CLR_file_truth = 'ADT_QC_CLR_count_GSE128639_20230115.csv'
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
#protein_list = pd.read_csv('union_50_proteins_all_4_external_val_sets.csv')['x'].tolist()
training_RNA_combined = pd.read_csv('all_top20_RNA_feature_combined_input_RNA_ref_combined_6_datasets_all_training_samples_20230115.csv', index_col = 'Unnamed: 0')
RNA_test = pd.read_csv(data_denoised_estimate_QC_file, index_col = 'Unnamed: 0')
test_ADT_truth = pd.read_csv(ADT_CLR_file_truth, index_col = 'Unnamed: 0')
protein_list = test_ADT_truth.index.tolist()

#Training:
param={}
for p in protein_list:
    RNA_protein_corr = pd.read_csv(p + '_RNA_feature_selection_cor_ref_combined_6_datasets_all_training_samples_20230115.csv', index_col = 'Unnamed: 0')
    shared_RNA_features = list(set(RNA_protein_corr.index) & set(RNA_test.index))
    RNA_protein_corr = RNA_protein_corr.loc[shared_RNA_features, :]
    rna_sort = RNA_protein_corr.sort_values('corr', ascending=False)
    high_corr_mrna = rna_sort[0:20].index.tolist()
    train_cell_id = pd.read_csv(p + '_cell_features_training_combined_6_sets_all_training_samples_cell_id_20230115.csv')['train_cell_id'] #Import training set cell id
    x = training_RNA_combined.loc[high_corr_mrna, train_cell_id].T
    y=  pd.read_csv(p + '_ADT_CLR_truth_ref_combined_6_datasets_all_training_samples_20230115.csv', index_col = 'Unnamed: 0')[train_cell_id].T
    #train RandomForest regressor:
    rf2 = RandomForestRegressor(**param)
    rf2.fit(x, y)
    with open(p + '_' + dataset_id, 'wb') as f:
        pickle.dump(rf2, f)

