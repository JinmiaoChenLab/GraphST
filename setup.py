from setuptools import Command, find_packages, setup

__lib_name__ = "GraphST"
__lib_version__ = "1.1.1"
__description__ = "Spatially informed clustering, integration, and deconvolution of spatial transcriptomics with GraphST"
__url__ = "https://github.com/JinmiaoChenLab/GraphST"
__author__ = "Yahui Long"
__author_email__ = "longyh@immunol.a-star.edu.sg"
__license__ = "MIT"
__keywords__ = ["Spatial transcriptomics", "Graph self-supervised contrastive learning", "Spatial clustering", "Cell type deconvolution", "Batch correction", "Single-cell RNA-seq"]
__requires__ = ["requests",]

with open("README.rst", "r", encoding="utf-8") as f:
    __long_description__ = f.read()

setup(
    name = __lib_name__,
    version = __lib_version__,
    description = __description__,
    url = __url__,
    author = __author__,
    author_email = __author_email__,
    license = __license__,
    packages = ["GraphST"],
    install_requires = __requires__,
    zip_safe = False,
    include_package_data = True,
    long_description = """GraphST is a versatile graph self-supervised contrastive learning model that incorporates spatial location information and gene expression profiles to accomplish three key tasks, spatial clustering, spatial transcriptomics (ST) data integration, and single-cell RNA-seq (scRNA-seq) transfer onto ST. GraphST combines graph neural networks (GNNs) with self-supervised contrastive learning to learn spot representations in the ST data by modeling gene expressions and spatial locaiton information. After the representation learning, the non-spatial alignment algorithm is used to cluster the spots into different spatial domains. Each cluster is regarded as a spatial domain, containing spots with similar gene expression profiles and spatially proximate. GraphST can jointly analyze multiple ST samples while correcting batch effects, which is achieved by smoothing features between spatially adjacent spots across samples. For the scRNA-seq transfer onto ST data, a mapping matrix is trained via an augmentation-free contrastive learning mechanism, where the similarity of spatially adjacent spots are maximized while those of spatially non-adjacent spots are minimized. With the learned mapping matrix, arbitrary cell attributes (e.g., cell type and sample type) can be flexibly projected onto spatial space. """,
    long_description_content_type="text/markdown"
)
