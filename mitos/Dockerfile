FROM continuumio/miniconda:4.7.12
MAINTAINER <philipp.resl@uni-graz.at>

RUN conda install -y -c bioconda mitos=1.0.5
RUN conda clean -a -y 
WORKDIR /data

CMD ["runmitos.py"]
