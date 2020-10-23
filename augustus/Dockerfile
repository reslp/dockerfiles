FROM continuumio/miniconda:4.7.10
MAINTAINER <philipp.resl@uni-graz.at>

ENV PATH /opt/conda/bin:$PATH

RUN conda config --append channels bioconda
RUN conda config --append channels conda-forge
RUN conda config --append channels anaconda
RUN conda install -c bioconda augustus=3.3.3
RUN conda clean -a -y

WORKDIR /data

CMD ["augustus", "-h"]
