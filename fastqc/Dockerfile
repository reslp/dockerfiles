FROM continuumio/miniconda:4.7.10

MAINTAINER <resl.philipp@bio.lmu.de>

RUN conda config --set notify_outdated_conda false
RUN conda install -y -c bioconda fastqc=0.11.7
RUN conda clean -a -y

CMD ["fastqc"]

