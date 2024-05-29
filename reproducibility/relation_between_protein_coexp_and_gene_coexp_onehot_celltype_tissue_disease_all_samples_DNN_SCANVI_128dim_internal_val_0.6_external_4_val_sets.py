# source /opt/rh/rh-python36/enable
# source /home/stat/zhouzilu/project/zhouzilu/expressionGAN/.env/bin/activate
import numpy as np
import torch.nn as nn
import pandas as pd
import matplotlib.pyplot as plt
import torch.nn.functional as F
import torch
import torch.optim as optim
import sys
from scipy import stats
import os
import random
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import Normalizer
from sklearn.preprocessing import normalize
import copy
import time

from sklearn.linear_model import LinearRegression
import pickle
from numpy.linalg import norm
from sklearn.metrics.pairwise import cosine_similarity



#--------------------------------------------------------------------------------------
#changeable parameters:

#GSE128639:
use_n_ensemble = 8

threshold_DNN_internal_val_acc = '0.6'
dataset_id = 'GSE128639'
RNA_SCANVI_latent_representation_file = 'GSE128639_query_RNA_latent_representations_SCANVI_128dim_impute_HVG_zeros_trained_6_datasets_20230115.csv'

cell_type_file = 'celltype_cell_id_external_val_set_GSE128639_SingleR_20230115.csv'
tissue_type = 'bone marrow'
disease_type = 'healthy'


training_gene_coexp_matrix_file = 'all_celltypes_combined_6_training_sets_pearson_gene_cor_mat.csv'
test_gene_coexp_matrix_file = 'gene_cor_mat_pearson_GSE128639_.csv'
ADT_QC_CLR_norm_file = 'ADT_QC_CLR_count_GSE128639_20230115.csv'
test_set_protein_gene_name_file = 'all_25_protein_gene_name_GSE128639_BM_20230115.csv'

save_gix_test_set_file_path = '/few_shot_learning_protein_specific_onehot_celltype_tissue_disease_DNN_ensemble/'
save_dot_product_protein_features_path = '/zero_shot_learning_protein_specific_onehot_celltype_tissue_disease_DNN_ensemble/'
save_linear_model_pred_vs_truth_y_unseen_file_path = '/zero_shot_learning_protein_specific_onehot_celltype_tissue_disease_DNN_ensemble/'

#--------------------------------------------------------------------------------------
#changeable parameters:

#healthy pancreas GSM5025059:
use_n_ensemble = 8

threshold_DNN_internal_val_acc = '0.6'
dataset_id = 'GSM5025059'
RNA_SCANVI_latent_representation_file = 'pancreas_normal_GSM5025059_query_RNA_latent_representations_SCANVI_latent_128dim_impute_HVG_zeros_trained_6_datasets_20230115.csv'

cell_type_file = 'celltype_cell_id_external_val_set_GSM5025059_SingleR_20230115.csv'
tissue_type = 'pancreas'
disease_type = 'healthy'

training_gene_coexp_matrix_file = 'all_celltypes_combined_6_training_sets_pearson_gene_cor_mat.csv'
test_gene_coexp_matrix_file = 'gene_cor_mat_pearson_GSM5025059_pancreas_normal__.csv'
ADT_QC_CLR_norm_file = 'ADT_QC_CLR_count_pancreas_normal_GSM5025059_20230115.csv'
test_set_protein_gene_name_file = 'all_13_protein_gene_name_GSM5025059_20230115.csv'

save_gix_test_set_file_path = '/few_shot_learning_protein_specific_onehot_celltype_tissue_disease_DNN_ensemble/'
save_dot_product_protein_features_path = '/zero_shot_learning_protein_specific_onehot_celltype_tissue_disease_DNN_ensemble/'
save_linear_model_pred_vs_truth_y_unseen_file_path = '/zero_shot_learning_protein_specific_onehot_celltype_tissue_disease_DNN_ensemble/'

#--------------------------------------------------------------------------------------
#changeable parameters:

#pancreas_pancreatitis_GSM5025052:
use_n_ensemble = 8

threshold_DNN_internal_val_acc = '0.6'
dataset_id = 'GSM5025052'
RNA_SCANVI_latent_representation_file = 'pancreas_pancreatitis_GSM5025052_query_RNA_latent_representations_SCANVI_latent_128dim_impute_HVG_zeros_trained_6_datasets_20230115.csv'

cell_type_file = 'celltype_cell_id_external_val_set_GSM5025052_SingleR_20230115.csv'
tissue_type = 'pancreas'
disease_type = 'pancreatitis'

training_gene_coexp_matrix_file = 'all_celltypes_combined_6_training_sets_pearson_gene_cor_mat.csv'
test_gene_coexp_matrix_file = 'gene_cor_mat_pearson_GSM5025052_pancreas_pancreatitis__.csv'
ADT_QC_CLR_norm_file = 'ADT_QC_CLR_count_pancreas_pancreatitis_GSM5025052_20230115.csv'
test_set_protein_gene_name_file = 'all_13_protein_gene_name_GSM5025052_20230115.csv'

save_gix_test_set_file_path = '/few_shot_learning_protein_specific_onehot_celltype_tissue_disease_DNN_ensemble/'
save_dot_product_protein_features_path = '/zero_shot_learning_protein_specific_onehot_celltype_tissue_disease_DNN_ensemble/'
save_linear_model_pred_vs_truth_y_unseen_file_path = '/zero_shot_learning_protein_specific_onehot_celltype_tissue_disease_DNN_ensemble/'

#--------------------------------------------------------------------------------------
#changeable parameters:

#COVID BALF GSM5093918:
use_n_ensemble = 8

threshold_DNN_internal_val_acc = '0.6'
dataset_id = 'GSM5093918'
RNA_SCANVI_latent_representation_file = 'COVID_lung_GSM5093918_Sar01L_CD3pos_RNA_latent_representations_SCANVI_latent_128dim_impute_HVG_zeros_trained_6_datasets_20230115.csv'

cell_type_file = 'celltype_cell_id_external_val_set_GSM5093918_SingleR_20230115.csv'
tissue_type = 'BALF'
disease_type = 'COVID'

training_gene_coexp_matrix_file = 'all_celltypes_combined_6_training_sets_pearson_gene_cor_mat.csv'
test_gene_coexp_matrix_file = 'gene_cor_mat_pearson_GSE167118_lung_COVID__.csv'
ADT_QC_CLR_norm_file = 'ADT_QC_CLR_count_lung_COVID_GSM5093918_20230115.csv'
test_set_protein_gene_name_file = 'all_39_protein_gene_name_GSM5093918_20230115.csv'

save_gix_test_set_file_path = '/few_shot_learning_protein_specific_onehot_celltype_tissue_disease_DNN_ensemble/'
save_dot_product_protein_features_path = '/zero_shot_learning_protein_specific_onehot_celltype_tissue_disease_DNN_ensemble/'
save_linear_model_pred_vs_truth_y_unseen_file_path = '/zero_shot_learning_protein_specific_onehot_celltype_tissue_disease_DNN_ensemble/'

#--------------------------------------------------------------------------------------
#Prepare zero-shot learning:
#-------------------------------------------------------------
#Get each one of the ensemble members' ~10000-dim normalized original features (i.e., Z1, Z2, ... Z22 for P1, P2, ... P22, respectively), save them:
all_training_protein_list = pd.read_csv('protein_gene_names_selected_final_ensemble members_from_combined_6_training_sets_DNN_onehot_celltype_tissue_disease_SCANVI_128dim_internal_val_threshold_' + str(threshold_DNN_internal_val_acc) + '_20230115.csv')['consistent_protein_name']

x_training_set = pd.read_csv('all_celltypes_combined_6_training_sets_pearson_gene_cor_mat.csv', index_col = 'Unnamed: 0')

x = pd.read_csv(test_gene_coexp_matrix_file, index_col = 'Unnamed: 0')

all_trainable_proteins_gene_names_6_training_sets = pd.read_csv('protein_gene_names_union_289_DNNs_from_combined_6_training_sets_DNN_onehot_celltype_tissue_disease_SCANVI_128dim_internal_val_threshold_' + str(threshold_DNN_internal_val_acc) + '_20230115.csv')

all_trainable_proteins_gene_names_6_training_sets.index = all_trainable_proteins_gene_names_6_training_sets['consistent_protein_name']
all_trainable_proteins_gene_names_6_training_sets = all_trainable_proteins_gene_names_6_training_sets.loc[all_training_protein_list, :]

all_trainable_proteins_gene_names_6_training_sets_2 = all_trainable_proteins_gene_names_6_training_sets
all_trainable_proteins_gene_names_6_training_sets_2.index = all_trainable_proteins_gene_names_6_training_sets_2['gene_name']

#Select ensemble members with corresponding genes existing in external validation set's coexpression file:
shared_gene_ensemble_members_test_set_coexp = list(set(x.columns) & set(all_trainable_proteins_gene_names_6_training_sets_2['gene_name']))
use_all_training_protein_list = list(all_trainable_proteins_gene_names_6_training_sets_2.loc[shared_gene_ensemble_members_test_set_coexp, 'consistent_protein_name'])
#Select ensemble members with corresponding genes existing in external validation set's coexpression file

shared_training_test_set_coexp_features = list(set(x_training_set.columns) & set(x.columns))
z = x[shared_training_test_set_coexp_features]

z_ensemble = z.loc[shared_gene_ensemble_members_test_set_coexp, :]
z_ensemble.index = use_all_training_protein_list


#-------------------------------------------------------------
#Get each one of the tested unseen proteins' ~10000-dim normalized original features (i.e., Z1, Z2, ... Z22 for P1, P2, ... P22, respectively), save them:

all_test_proteins_gene_names = pd.read_csv(test_set_protein_gene_name_file, index_col = 'Unnamed: 0') 
#all_test_proteins_gene_names.index = all_test_proteins_gene_names['protein_names']
all_test_proteins_gene_names.index = all_test_proteins_gene_names.loc[:, 'protein_name']

all_test_proteins_gene_names_2 = all_test_proteins_gene_names
all_test_proteins_gene_names_2.index = all_test_proteins_gene_names.loc[:, 'gene_name']


shared_gene_unseen_proteins_test_set_coexp = list(set(x.columns) & set(all_test_proteins_gene_names.iloc[:, 1]))
use_all_protein_list = all_test_proteins_gene_names_2.loc[shared_gene_unseen_proteins_test_set_coexp, :]
use_all_protein_list.index = use_all_protein_list.iloc[:, 0]
use_all_protein_list = use_all_protein_list.iloc[:, 0] #unseen proteins

if dataset_id == 'GSE128639':
    use_all_protein_list = use_all_protein_list.drop({'CD45RO'})

if dataset_id == 'GSM5025059':
    use_all_protein_list = use_all_protein_list.drop({'CD45'})

if dataset_id == 'GSM5025052':
    use_all_protein_list = use_all_protein_list.drop({'CD45'})

if dataset_id == 'GSM5093918':
    use_all_protein_list = use_all_protein_list.drop({'CD45', 'CD45RO'})

z_unseen_proteins = z.loc[shared_gene_unseen_proteins_test_set_coexp, :]
z_unseen_proteins.index = use_all_protein_list

use_common_genes = list(set(z_unseen_proteins.columns) & set(z_ensemble.columns))
#df_cosine_similarity = cosine_similarity(z_ensemble[use_common_genes], z_unseen_proteins[use_common_genes])
#df_cosine_similarity = pd.DataFrame(df_cosine_similarity, index = z_ensemble.index, columns = z_unseen_proteins.index)

df_gene_coexp_cosine_similarity = cosine_similarity(z_unseen_proteins[use_common_genes], z_unseen_proteins[use_common_genes])
df_gene_coexp_cosine_similarity = pd.DataFrame(df_gene_coexp_cosine_similarity, index = z_unseen_proteins.index, columns = z_unseen_proteins.index)

ADT2 = pd.read_csv(ADT_QC_CLR_norm_file, index_col = 'Unnamed: 0')
all_protein_list = use_all_protein_list
ADT = ADT2   
ADT = pd.DataFrame(ADT.loc[all_protein_list, :])
df_protein_exp_cosine_similarity = cosine_similarity(ADT, ADT)
df_protein_exp_cosine_similarity = pd.DataFrame(df_protein_exp_cosine_similarity, index = ADT.index, columns = ADT.index)


df = pd.DataFrame({'df_protein_exp_cosine_similarity':df_protein_exp_cosine_similarity.values.ravel().tolist(), 'df_gene_coexp_cosine_similarity': df_gene_coexp_cosine_similarity.values.ravel().tolist()})
df.corr()


def softmax(x):
    """Compute softmax values for each sets of scores in x."""
    return np.exp(x) / np.sum(np.exp(x), axis=1)


softmax_df_gene_coexp_cosine_similarity = softmax(df_gene_coexp_cosine_similarity)
softmax_df_protein_exp_cosine_similarity = softmax(df_protein_exp_cosine_similarity)


df3 = pd.DataFrame({'df_protein_exp_cosine_similarity':softmax_df_protein_exp_cosine_similarity.values.ravel().tolist(), 'df_gene_coexp_cosine_similarity': softmax_df_gene_coexp_cosine_similarity.values.ravel().tolist()})

#GSE128639:
#softmax_df_gene_coexp_cosine_similarity.to_csv(dataset_id + '_softmax_df_gene_coexp_cosine_similarity_20240504.csv')
#softmax_df_protein_exp_cosine_similarity.to_csv(dataset_id + '_softmax_df_protein_exp_cosine_similarity_20240504.csv')
#healthy pancreas GSM5025059:
#softmax_df_gene_coexp_cosine_similarity.to_csv(dataset_id + '_softmax_df_gene_coexp_cosine_similarity_20240504.csv')
#softmax_df_protein_exp_cosine_similarity.to_csv(dataset_id + '_softmax_df_protein_exp_cosine_similarity_20240504.csv')
#pancreas_pancreatitis_GSM5025052:
#softmax_df_gene_coexp_cosine_similarity.to_csv(dataset_id + '_softmax_df_gene_coexp_cosine_similarity_20240504.csv')
#softmax_df_protein_exp_cosine_similarity.to_csv(dataset_id + '_softmax_df_protein_exp_cosine_similarity_20240504.csv')
#COVID BALF GSM5093918:
softmax_df_gene_coexp_cosine_similarity.to_csv(dataset_id + '_softmax_df_gene_coexp_cosine_similarity_20240504.csv')
softmax_df_protein_exp_cosine_similarity.to_csv(dataset_id + '_softmax_df_protein_exp_cosine_similarity_20240504.csv')









        

