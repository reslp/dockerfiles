FROM continuumio/miniconda3:4.7.10
MAINTAINER <resl.philipp@bio.lmu.de>

ENV PATH /opt/conda/bin:$PATH

RUN conda config --append channels bioconda && \
	conda config --append channels conda-forge && \
	conda config --append channels anaconda && \
	conda install -c bioconda pandas=0.25.1 && \
	conda install -c bioconda concoct=1.1 && \
	conda clean -a -y

CMD ["concoct", "-h"]