FROM reslp/mamba:0.15.3
MAINTAINER <philipp.resl@uni-graz.at>

RUN mamba install -y -c bioconda colord=1.1.0

WORKDIR /data

CMD ["colord"]
