FROM reslp/mamba:0.15.3
MAINTAINER <philipp.resl@uni-graz.at>

RUN mamba install -y -c bioconda bmge=1.12

WORKDIR /DATA

CMD ["bmge"]
