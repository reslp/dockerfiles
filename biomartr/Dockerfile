FROM continuumio/miniconda:4.7.10
MAINTAINER <resl.philipp@bio.lmu.de>

ENV PATH /opt/conda/bin:$PATH
RUN apt-get update && apt-get install -y build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev curl
RUN conda config --append channels bioconda && \
	conda config --append channels conda-forge && \
	conda config --append channels anaconda && \
	conda install -c bioconda bioconductor-biomart
RUN conda install -c conda-forge r-devtools
RUN conda install -c bioconda bioconductor-biostrings
RUN conda clean -a -y
RUN R -e "install.packages('biomartr', dependencies = TRUE, repos='http://cran.rstudio.com/')"

CMD ["/bin/bash"]