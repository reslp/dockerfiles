FROM continuumio/miniconda3:4.7.10
MAINTAINER <resl.philipp@bio.lmu.de>

ENV PATH /opt/conda/bin:$PATH

RUN conda config --append channels bioconda && \
	conda config --append channels conda-forge && \
	conda config --append channels anaconda && \
	conda config --set notify_outdated_conda false && \
	conda install -c anaconda python=3.6.8 && \
	conda install -c bioconda multiqc=1.6 && \
	conda clean -a -y

ENTRYPOINT ["multiqc"]
CMD ["."]