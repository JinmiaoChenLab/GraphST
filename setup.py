from setuptools import Command, find_packages, setup

with open("README.rst", "r", encoding="utf-8") as f:
    __description__ = f.read()

setup(
    name = "DeepST",
    version = "1.0.0",
    description = "DeepST: A versatile graph contrastive learning framework for spatially informed clustering, integration and transfer of spatial transcriptomics",
    url = "https://github.com/longyahui/DeepST",
    author = "Yahui Long",
    author_email = "longyh@immunol.a-star.edu.sg",
    license = "MIT",
    packages = ['DeepST'],
    install_requires = ["torch=1.8.0"],
    zip_safe = False,
    include_package_data = True,
    long_description = __description__
)
