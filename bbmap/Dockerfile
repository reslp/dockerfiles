FROM reslp/mamba:0.7.0
MAINTAINER <philipp.resl@uni-graz.at>

RUN mamba install -y -c bioconda bbmap=38.90

WORKDIR /data

CMD ["bbmap.sh"]
