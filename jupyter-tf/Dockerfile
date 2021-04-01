FROM jupyter/all-spark-notebook:latest
# this is simply extending the official jupyter/all-spark-notebook container with tensorflow 2

MAINTAINER <philipp.resl@uni-graz.at>

ENV PATH /opt/conda/bin:$PATH

RUN conda config --append channels bioconda
RUN conda config --append channels conda-forge
RUN conda config --append channels anaconda
RUN conda config --set notify_outdated_conda false

RUN conda install mamba
RUN mamba install tensorflow=2.4.1
