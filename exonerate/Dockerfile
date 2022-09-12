FROM reslp/mamba:0.21.1
MAINTAINER <philipp.resl@uni-graz.at>

RUN mamba install -c bioconda exonerate=2.4.0

WORKDIR /data

CMD ["exonerate"]
