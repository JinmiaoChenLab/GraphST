import numpy as np
import pandas as pd
from sklearn import metrics
from sklearn.decomposition import PCA


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

def clustering(adata, n_clusters=7, radius=50, key='emb', threshold=0.06, data_type='10X', refinement=True):
    """\
    Spatial clustering based the learned representation.

    Parameters
    ----------
    adata : anndata
        AnnData object of scanpy package.
    n_clusters : int, optional
        The number of clusters. The default is 7.
    radius : int, optional
        The number of neighbors considered during refinement. The default is 50.
    key : string, optional
        The key of the learned representation in adata.obsm. The default is 'emb'.
    threshold : float, optional
        Cutoff for selecting the final labels. For 10X Visium data, the model is trained twice,
        i.e., with and without penalty terms. As a result,  two clustering label results corresponding to the
        two-time training are generated. The final clustering label is determined by Silhouette score. 
        The default is 0.06.
    data_type : string, optional
        Data type of input spatial data. The current version of DeepST can be applied to different ST data,
        including 10X Visium, Stereo-seq, and Slide-seqV2.
    refinement : bool, optional
        Refine the predicted labels or not. The default is True.

    Returns
    -------
    None.

    """
    
    pca = PCA(n_components=20, random_state=42) 
    
    if data_type == '10X':
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
    
       # Silhouette
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
             
    elif data_type in ['Stereo', 'SlideV2']:
         embedding = pca.fit_transform(adata.obsm['emb'].copy())
         adata.obsm['emb_pca'] = embedding
         adata = mclust_R(adata, used_obsm='emb_pca', num_cluster=n_clusters)
         adata.obs['label'] = adata.obs['mclust'] 
       
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

def extract_top_value(map_matrix, percent = 0.1): 
    '''\
    Filter out cells with low mapping probability

    Parameters
    ----------
    map_matrix : array
        Mapped matrix with m spots and n cells.
    percent : float, optional
        The percentage of cells to retain. The default is 0.1.

    Returns
    -------
    output : array
        Filtered mapped matrix.

    '''

    #retain top 1% values for each spot
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
    mat = pd.DataFrame(zeros, index=adata_sc.obs_names, columns=cell_type)
    for cell in list(adata_sc.obs_names):
        ctype = adata_sc.obs.loc[cell, label]
        mat.loc[cell, str(ctype)] = 1
    #res = mat.sum()
    return mat

def project_cell_to_spot(adata, adata_sc):
    '''
    Project cell types onto ST data using mapped matrix in adata.obsm

    Parameters
    ----------
    adata : anndata
        AnnData object of spatial data.
    adata_sc : anndata
        AnnData object of scRNA-seq reference data.

    Returns
    -------
    None.

    '''
    
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
