FROM reslp/mamba:0.7.0
MAINTAINER <philipp.resl@uni-graz.at>

RUN mamba install purge_haplotigs=1.1.1

WORKDIR /data 

CMD ["purge_haplotigs"]
