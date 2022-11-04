"""
Setup
-------------------------------------------------------------------------------
The best way to install NanoDiP is to create a virtual environment and then
install via pip, after updating setuptools i.e.:

virtualenv -p {path_to_python} {virtual_env_name}
source {virtual_env_path}/bin/activate
pip install -U setuptools
pip install -e .

"""

import sys

from setuptools import find_packages
from setuptools import setup

sys.path.insert(0, "/applications/nanodip")
from nanodip.config import __version__ as VERSION

# Load the README file.
with open(file="README.md", mode="r") as f:
    long_description = f.read()

setup(
    name="nanodip",
    author="J. Hench, C. Hultschig, J. Brugger, and S. Frank, Neuropathology IfP Basel",
    version=VERSION,
    description="NanoDiP analyzes Nanopore sequencing and methylation data for tumor classification.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/brj0/nanodip",
    install_requires=[
        "cherrypy",
        # The current version of grpcio (1.50.0) interferes with minknow_api
        # and causes subprocess.Popen (used in nanodip) to freeze occasionally.
        "grpcio==1.44.0",
        "ipython",
        "jinja2",
        "kaleido",
        # CAUTION: Requires a patched version of minknow api, file
        # `[VENV]/lib/python3.7/site-packages/minknow_api/tools/protocols.py`.
        # Without the patch, the generated fast5 sequencing data will be
        # unreadable with f5c or nanopolish (wrong compression algorithm, which
        # is the default in the MinKNOW backend).
        "minknow_api@git+https://github.com/neuropathbasel/minknow_api",
        "nbformat",
        "numpy",
        "openpyxl",
        "pandas",
        "pdf2image",
        "plotly",
        "psutil",
        "pysam",
        "reportlab==3.6.1",
        "scikit_learn",
        "scipy",
        "statsmodels",
        "tqdm",
        "umap-learn",
        "xhtml2pdf==0.2.5",
    ],
    extras_require={
        "cupy": ["cupy"],
    },
    packages=find_packages(),
    # Include package data such as images.
    include_package_data=True,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: GPU :: NVIDIA CUDA",
        "Framework :: CherryPy",
        "Framework :: Jupyter :: JupyterLab",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
    ],
    entry_points={"console_scripts": ["nanodip=nanodip.main:start_nanodip"]},
)
