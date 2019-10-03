FROM continuumio/miniconda3:4.7.10
MAINTAINER <resl.philipp@bio.lmu.de>

ENV PATH /opt/conda/bin:$PATH
RUN conda config --set notify_outdated_conda false
RUN conda config --append channels bioconda
RUN conda config --append channels conda-forge
RUN conda config --append channels anaconda
RUN conda install -c conda-forge biopython
ADD extract_contigs_from_blobtools.py /usr/bin/extract_contigs_from_blobtools.py
RUN chmod a+x /usr/bin/extract_contigs_from_blobtools.py
WORKDIR /data/

ENTRYPOINT ["/usr/bin/extract_contigs_from_blobtools.py"]