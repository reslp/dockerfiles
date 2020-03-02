FROM continuumio/miniconda:4.7.10
MAINTAINER <philipp.resl@uni-graz.at>

ENV PATH /opt/conda/bin:$PATH

RUN conda config --append channels bioconda && \
	conda config --append channels conda-forge && \
	conda config --append channels anaconda && \
	conda install -c bioconda binsanity=0.3.4 && \
	conda clean -a -y

CMD ["Binsanity-profile"]