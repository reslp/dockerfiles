FROM continuumio/miniconda:4.7.10
MAINTAINER <resl.philipp@bio.lmu.de>

ENV PATH /opt/conda/bin:$PATH

RUN conda config --append channels bioconda
RUN conda config --append channels conda-forge
RUN conda config --append channels anaconda
RUN conda install -c bioconda bowtie2=2.3.5
RUN conda clean -a -y

CMD ["bowtie2", "-h"]