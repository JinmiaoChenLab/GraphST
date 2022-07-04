from sklearn.decomposition import PCA
import pandas as pd
import seaborn as sns
import os
import numpy as np
from sklearn import metrics
import matplotlib.pyplot as plt
import matplotlib as mpl
import scanpy as sc

os.environ['R_HOME'] = '/scbio4/tools/R/R-4.0.3_openblas/R-4.0.3'

def mclust_R(adata, num_cluster, modelNames='EEE', used_obsm='emb_pca', random_seed=2020):
    """\
    Clustering using the mclust algorithm.
    The parameters are the same as those in the R package mclust.
    """
    
    np.random.seed(random_seed)
    import rpy2.robjects as robjects
    robjects.r.library("mclust")

    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    r_random_seed = robjects.r['set.seed']
    r_random_seed(random_seed)
    rmclust = robjects.r['Mclust']
    
    res = rmclust(rpy2.robjects.numpy2ri.numpy2rpy(adata.obsm[used_obsm]), num_cluster, modelNames)
    mclust_res = np.array(res[-2])

    adata.obs['mclust'] = mclust_res
    adata.obs['mclust'] = adata.obs['mclust'].astype('int')
    adata.obs['mclust'] = adata.obs['mclust'].astype('category')
    return adata

def clustering(adata, n_clusters=7, radius=50, key='emb', threshold=0.06, refinement=True):
    pca = PCA(n_components=20, random_state=42) 
    
    # clustering 1
    embedding = pca.fit_transform(adata.obsm['emb'].copy())
    adata.obsm['emb_pca'] = embedding
    adata = mclust_R(adata, used_obsm='emb_pca', num_cluster=n_clusters)
    adata.obs['label'] = adata.obs['mclust']
    new_type = refine_label(adata, radius, key='label')
    adata.obs['label_refined'] = new_type
    
    # clustering 2
    embedding = pca.fit_transform(adata.obsm['emb_reg'].copy())
    adata.obsm['emb_reg_pca'] = embedding
    adata = mclust_R(adata, used_obsm='emb_reg_pca', num_cluster=n_clusters)
    adata.obs['label_reg'] = adata.obs['mclust']
    new_type = refine_label(adata, radius, key='label_reg')
    adata.obs['label_reg_refined'] = new_type
    
    # Calculate Silhouette score
    SIL = metrics.silhouette_score(adata.obsm['emb_pca'], adata.obs['label'], metric='euclidean')
    SIL_reg = metrics.silhouette_score(adata.obsm['emb_reg_pca'], adata.obs['label_reg'], metric='euclidean')
    
    if abs(SIL-SIL_reg) > threshold and SIL_reg > SIL:
       if refinement: 
          adata.obs['domain'] = adata.obs['label_reg_refined']
       else:   
          adata.obs['domain'] = adata.obs['label_reg']
    else:
       if refinement: 
          adata.obs['domain'] = adata.obs['label_refined']
       else:
          adata.obs['domain'] = adata.obs['label'] 
       
def refine_label(adata, radius=50, key='label'):
    n_neigh = radius
    new_type = []
    old_type = adata.obs[key].values
    
    #read distance
    if 'distance_matrix' not in adata.obsm.keys():
        raise ValueError("Distance matrix is not existed!")
    distance = adata.obsm['distance_matrix'].copy()
           
    n_cell = distance.shape[0]
    
    for i in range(n_cell):
        vec  = distance[i, :]
        index = vec.argsort()
        neigh_type = []
        for j in range(1, n_neigh+1):
            neigh_type.append(old_type[index[j]])
        max_type = max(neigh_type, key=neigh_type.count)
        new_type.append(max_type)
        
    new_type = [str(i) for i in list(new_type)]    
    #adata.obs['label_refined'] = np.array(new_type)
    
    return new_type

def extract_top_value(map_matrix, percent = 0.1): # 0.05
    """
    map_matrix: projection matrix with m cells and n spots.
    """
    #retain top 5% values for each cell
    top_k  = percent * map_matrix.shape[1]
    output = map_matrix * (np.argsort(np.argsort(map_matrix)) >= map_matrix.shape[1] - top_k)
    
    return output 

def construct_cell_type_matrix(adata_sc):
    label = 'cell_type'
    n_type = len(list(adata_sc.obs[label].unique()))
    zeros = np.zeros([adata_sc.n_obs, n_type])
    cell_type = list(adata_sc.obs[label].unique())
    cell_type = [str(s) for s in cell_type]
    cell_type.sort()
    print('cell_type:', cell_type)
    mat = pd.DataFrame(zeros, index=adata_sc.obs_names, columns=cell_type)
    for cell in list(adata_sc.obs_names):
        ctype = adata_sc.obs.loc[cell, label]
        mat.loc[cell, str(ctype)] = 1
    print('shape of mat:', mat.shape)
    res = mat.sum()
    print(res)
    return mat

def project_cell_to_spot(adata, adata_sc):
    plt.rcParams['font.sans-serif'] = ['Times New Roman']
    
    # read map matrix 
    map_matrix = adata.obsm['map_matrix']   # spot x cell
   
    # extract top-k values for each spot
    map_matrix = extract_top_value(map_matrix) # filtering by spot
    
    # construct cell type matrix
    matrix_cell_type = construct_cell_type_matrix(adata_sc)
    matrix_cell_type = matrix_cell_type.values
       
    # projection by spot-level
    matrix_projection = map_matrix.dot(matrix_cell_type)
   
    # rename cell types
    cell_type = list(adata_sc.obs['cell_type'].unique())
    cell_type = [str(s) for s in cell_type]
    cell_type.sort()
    #cell_type = [s.replace(' ', '_') for s in cell_type]
    df_projection = pd.DataFrame(matrix_projection, index=adata.obs_names, columns=cell_type)  # spot x cell type
    
    #normalize by row (spot)
    df_projection = df_projection.div(df_projection.sum(axis=1), axis=0).fillna(0)

    #add projection results to adata
    adata.obs[df_projection.columns] = df_projection
    print('adata:', adata)
    
    with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
            
         ax = sc.pl.spatial(adata, cmap='magma',
                  color = cell_type,
                  #color = ['adjacent normal', 'solid tumor'],
                  ncols=5, size=1.0,
                  wspace = 0.1, hspace = 0.2,
                  #wspace = 0.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2',
                  show=True
                 )

          
