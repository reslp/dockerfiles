FROM continuumio/miniconda:4.7.12
MAINTAINER <philipp.resl@uni-graz.at>

RUN conda install -y -c bioconda mitos=2.0.8
RUN conda clean -a -y 
WORKDIR /data

CMD ["runmitos.py"]
