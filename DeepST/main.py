import scanpy as sc
from DeepST import DeepST
from utils import clustering
import pandas as pd
from sklearn import metrics
import torch
import multiprocessing as mp

#device = torch.device('cuda:3' if torch.cuda.is_available() else 'cpu')

n_clusters = 7
radius = 50
dataset = '151507'
print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   {}   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'.format(dataset))     
# read data
file_fold = '/home/yahui/Yahui/Projects/data/' + str(dataset)
adata = sc.read_visium(file_fold, count_file='filtered_feature_bc_matrix.h5', load_images=True)
adata.var_names_make_unique()

model = DeepST(adata)
adata = model.train_DeepST()

print('adata:', adata)
  
clustering(adata, n_clusters, radius, refinement=True)

# add ground_truth
df_meta = pd.read_csv(file_fold + '/metadata.tsv', sep='\t')
df_meta_layer = df_meta['layer_guess']
#df_meta_layer = df_meta['ground_truth']
adata.obs['ground_truth'] = df_meta_layer.values

# filter out NA nodes
adata = adata[~pd.isnull(adata.obs['ground_truth'])] 
        
# calculate ARI
ARI = metrics.adjusted_rand_score(adata.obs['domain'], adata.obs['ground_truth'])
#ARI_refined = metrics.adjusted_rand_score(adata.obs['domain_refined'], adata.obs['ground_truth'])
adata.uns['ARI'] = ARI
    
print('Dataset:', dataset)
print('ARI:', ARI)

# plotting clustering result
sc.pl.spatial(adata, img_key="hires", color=["domain"], title=['DeepST (ARI=%.4f)'%ARI], show=True)







