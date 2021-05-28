FROM reslp/mamba:0.7.0
MAINTAINER <philipp.resl@uni-graz.at>

RUN mamba install -y -c bioconda sra-tools=2.10.9
RUN printf '/LIBS/GUID = "%s"\n' `uuidgen` > /root/.ncbi/user-settings.mkfg

WORKDIR /data

CMD ["/bin/bash"]
