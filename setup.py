from setuptools import Command, find_packages, setup

__lib_name__ = "DeepST"
__lib_version__ = "1.0.0"
__description__ = "DeepST: A versatile graph contrastive learning framework for spatially informed clustering, integration and transfer of spatial transcriptomics"
__url__ = "https://github.com/longyahui/DeepST"
__author__ = "Yahui Long"
__author_email__ = "longyh@immunol.a-star.edu.sg"
__license__ = "MIT"
__keywords__ = ["spatial transcriptomics", "Deep learning", "Graph neural networks", "Contrastive learning", "Self-supervised learning"]
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
    packages = ['DeepST'],
    install_requires = __requires__,
    zip_safe = False,
    include_package_data = True,
    long_description = __long_description__
)
