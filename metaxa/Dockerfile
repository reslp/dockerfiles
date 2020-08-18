FROM continuumio/miniconda3:4.7.10

MAINTAINER <philipp.resl@uni-graz.at>

RUN conda config --append channels bioconda
RUN conda config --append channels conda-forge
RUN conda config --append channels anaconda

RUN conda config --set notify_outdated_conda false
RUN conda install -c bioconda metaxa=2.2
RUN conda install -c bioconda blast-legacy
RUN conda clean -a -y

RUN mkdir data
WORKDIR data

CMD ["metaxa2", "-h"]


