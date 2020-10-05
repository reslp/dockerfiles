FROM reslp/mamba:0.5.3
MAINTAINER <philipp.resl@uni-graz.at>

RUN mamba install -y -c bioconda abyss=2.2.5

WORKDIR /data

CMD ["abyss-pe", "-h"]
