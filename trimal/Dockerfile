FROM reslp/mamba:0.5.3
MAINTAINER <philipp.resl@uni-graz.at>

RUN mamba install -y -c bioconda trimal=1.4.1

WORKDIR /data

CMD ["trimal", "-h"]
