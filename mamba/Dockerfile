FROM continuumio/miniconda3:4.8.2
MAINTAINER <philipp.resl@uni-graz.at>

ENV PATH /opt/conda/bin:$PATH

RUN conda config --append channels bioconda && \
	conda config --append channels conda-forge && \
	conda config --append channels anaconda && \
	conda config --set notify_outdated_conda false && \
	conda install -y -c conda-forge mamba=0.7.7 && \
	conda clean -a -y

CMD ["mamba"]
