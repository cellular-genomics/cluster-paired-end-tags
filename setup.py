"""Package setup."""

import setuptools

version = "1.0.0"

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cluster_PETs",
    version=version,
    author="Jakub Lipiński",
    author_email="jakub@cellular-genomics.com",
    description="Python library implementing paired end tags (PETs) clustering algorithm",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/cellular-genomics/cluster-paired-end-tags",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    python_requires='>=3.6',
    install_requires=[
        "pandas>=1.0.5",
        "numba>=0.50.1"
    ],
    entry_points={
        'console_scripts': ['cluster_PETs=cluster_pets.cluster_PETs:main'],
    }
)
