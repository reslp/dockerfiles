FROM reslp/mamba:0.7.0
MAINTAINER <philipp.resl@uni-graz.at>

RUN mamba install -y -c bioconda quicktree=2.5

WORKDIR /data

CMD ["quicktree"]
