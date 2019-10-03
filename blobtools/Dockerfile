FROM continuumio/miniconda3:4.7.10
MAINTAINER <resl.philipp@bio.lmu.de>


ENV PATH /opt/conda/bin:$PATH

RUN conda config --set notify_outdated_conda false
RUN conda config --append channels bioconda
RUN conda config --append channels conda-forge
RUN conda config --append channels anaconda
RUN conda install -c anaconda matplotlib docopt tqdm wget pyyaml git
RUN conda install -c bioconda pysam --update-deps
RUN conda clean -a -y
RUN git clone https://github.com/DRL/blobtools.git
WORKDIR /blobtools

WORKDIR /opt/db/
RUN wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
RUN tar zxf taxdump.tar.gz

WORKDIR /blobtools
RUN ./blobtools nodesdb --nodes /opt/db/nodes.dmp --names /opt/db/names.dmp


CMD ["-h"]
ENTRYPOINT ["./blobtools"]