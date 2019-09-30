FROM continuumio/miniconda:4.7.10

MAINTAINER <resl.philipp@bio.lmu.de>

ENV PATH /opt/conda/bin:$PATH

RUN conda config --append channels bioconda
RUN conda config --append channels conda-forge
RUN conda config --append channels anaconda

RUN conda config --set notify_outdated_conda false
RUN conda install -y -c bioconda busco=3.0.2
RUN conda clean -a -y

CMD ["run_busco", "-h"]