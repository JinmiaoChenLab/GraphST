from setuptools import Command, find_packages, setup

__lib_name__ = "GraphST"
__lib_version__ = "1.0.0"
__description__ = "Spatially informed clustering, integration, and deconvolution of spatial transcriptomics with GraphST"
__url__ = "https://github.com/longyahui/GraphST"
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
    packages = ['GraphST'],
    install_requires = __requires__,
    zip_safe = False,
    include_package_data = True,
    long_description = __long_description__
)
