FROM continuumio/miniconda:4.7.10
MAINTAINER <resl.philipp@bio.lmu.de>

ENV PATH /opt/conda/bin:$PATH

RUN conda config --append channels bioconda && \
	conda config --append channels conda-forge && \
	conda config --append channels anaconda && \
	conda install -c bioconda diamond=0.9.22 && \
	conda clean -a -y

CMD ["diamond", "version"]