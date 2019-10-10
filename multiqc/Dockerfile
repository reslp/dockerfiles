FROM continuumio/miniconda3:4.7.10

MAINTAINER <resl.philipp@bio.lmu.de>

ENV PATH /opt/conda/bin:$PATH

RUN conda config --append channels bioconda
RUN conda config --append channels conda-forge
RUN conda config --append channels anaconda

RUN conda config --set notify_outdated_conda false
RUN conda install -c anaconda python=3.6.8
RUN conda install -c bioconda multiqc=1.6
RUN conda clean -a -y

ENTRYPOINT ["multiqc"]
CMD ["."]