FROM continuumio/miniconda:4.7.10
MAINTAINER <resl.philipp@uni-graz.at>

ENV PATH /opt/conda/bin:$PATH

RUN conda config --append channels bioconda
RUN conda config --append channels conda-forge
RUN conda config --append channels anaconda
RUN conda install -c bioconda eggnog-mapper=2.0.1
RUN conda clean -a -y

CMD ["emapper.py"]