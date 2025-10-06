FROM reslp/mamba:0.21.1
MAINTAINER <philipp.resl@uni-graz.at>

RUN mamba install -y -c conda-forge magic-wormhole=0.12.0
WORKDIR /data

CMD ["wormhole"]

