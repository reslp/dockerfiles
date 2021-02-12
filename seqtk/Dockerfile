FROM reslp/mamba:0.7.0
MAINTAINER <philipp.resl@ui-graz.at>

RUN mamba install -c bioconda seqtk=1.3

WORKDIR /data

CMD ["seqtk"]
