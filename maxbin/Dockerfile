FROM continuumio/miniconda:4.7.10

MAINTAINER <resl.philipp@bio.lmu.de>

RUN conda install -c bioconda maxbin2=2.2.6
RUN conda clean -a -y

CMD ["run_MaxBin.pl"]