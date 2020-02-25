FROM continuumio/miniconda:4.7.10
MAINTAINER <resl.philipp@bio.lmu.de>

ENV PATH /opt/conda/bin:$PATH

RUN conda config --append channels bioconda && \
	conda config --append channels conda-forge && \
	conda config --append channels anaconda && \
	conda config --set notify_outdated_conda false && \
	conda install -y -c bioconda busco=3.0.2 && \
	conda clean -a -y
ENV AUGUSTUS_CONFIG_PATH /opt/conda/config
RUN chmod -R 777 /opt/conda/config
	
CMD ["run_busco", "-h"]