FROM reslp/mamba:0.7.0
MAINTAINER <philipp.resl@uni-graz.at>

RUN mamba install -y -c conda-forge -c bioconda edta=1.9.6 python=3.6 tensorflow=1.14 'h5py<3'

WORKDIR /data

CMD ["EDTA.pl"]
