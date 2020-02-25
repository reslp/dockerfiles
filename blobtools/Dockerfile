FROM continuumio/miniconda3:4.7.10
MAINTAINER <resl.philipp@bio.lmu.de>


ENV PATH /opt/conda/bin:$PATH

RUN conda config --set notify_outdated_conda false && \ 
	conda config --append channels bioconda && \
	conda config --append channels conda-forge  && \
	conda config --append channels anaconda  && \
	conda install -c anaconda matplotlib docopt tqdm wget pyyaml git  && \
	conda install -c bioconda pysam --update-deps  && \
	conda clean -a -y  && \
	wget https://github.com/DRL/blobtools/archive/blobtools_v1.1.1.tar.gz

RUN	tar zxf blobtools_v1.1.1.tar.gz && mv blobtools-blobtools_v1.1.1 blobtools

WORKDIR /opt/db/
RUN wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz && \
	tar zxf taxdump.tar.gz

WORKDIR /blobtools
RUN ./blobtools nodesdb --nodes /opt/db/nodes.dmp --names /opt/db/names.dmp

CMD ["-h"]
ENTRYPOINT ["./blobtools"]