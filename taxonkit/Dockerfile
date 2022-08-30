FROM reslp/mamba:0.21.1
MAINTAINER <philipp.resl@uni-graz.at>

RUN mamba install -c bioconda taxonkit=0.12.0

WORKDIR /data

CMD ["taxonkit"]
