FROM reslp/mamba:0.5.3
MAINTAINER <philipp.resl@uni-graz.at>

ENV PATH /opt/conda/bin:$PATH

RUN conda config --append channels bioconda && \
	conda config --append channels conda-forge && \
	conda config --append channels anaconda && \
	conda config --set notify_outdated_conda false && \
	mamba install -y -c bioconda busco=4.1.4 && \
	conda clean -a -y
ENV AUGUSTUS_CONFIG_PATH /opt/conda/config
RUN chmod -R 777 /opt/conda/config
	
CMD ["busco", "-h"]
