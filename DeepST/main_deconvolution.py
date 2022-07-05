import scanpy as sc
from DeepST import Train
from preprocess import preprocess, get_feature, construct_interaction, add_contrastive_label, fix_seed, filter_with_overlap_gene
from utils import project_cell_to_spot

#device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
device = 'cpu'

dataset = 'Human_Breast_Cancer'
random_seed = 50 
fix_seed(random_seed)

# read ST data
file_fold = '/home/yahui/Yahui/Projects/data/' + str(dataset)
adata = sc.read_visium(file_fold, count_file='filtered_feature_bc_matrix.h5', load_images=True)
adata.var_names_make_unique()

# preprocessing
preprocess(adata)

# built graph
construct_interaction(adata)
add_contrastive_label(adata)

# read scRNA daa
file_path = '/home/yahui/anaconda3/work/CellCluster_DEC/data/' + str(dataset) + '/scRNA_2.h5ad' 
adata_sc = sc.read(file_path)
adata_sc.var_names_make_unique()

# preprocessing
preprocess(adata_sc)

# find overlap genes
adata, adata_sc = filter_with_overlap_gene(adata, adata_sc)

# get features
get_feature(adata)

print('adata:', adata)
print('adata_sc:', adata_sc)

# Train model
model = Train(adata, adata_sc, epochs=1200, deconvolution=True, device=device)
adata, adata_sc = model.train_map()

# Project cells into spatial space
project_cell_to_spot(adata, adata_sc)




