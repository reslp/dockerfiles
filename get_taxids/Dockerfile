FROM continuumio/miniconda3:4.7.10
MAINTAINER <resl.philipp@bio.lmu.de>

ENV PATH /opt/conda/bin:$PATH

ADD get_taxids.py /usr/bin/get_taxids.py
RUN chmod a+x /usr/bin/get_taxids.py
#RUN wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid//prot.accession2taxid.gz
#RUN gunzip prot.accession2taxid.gz
WORKDIR /data/

CMD ["get_taxids.py"]