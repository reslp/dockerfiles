FROM reslp/mamba:0.7.0

MAINTAINER <philipp.resl@uni-graz.at>

RUN mamba install -c bioconda orthofinder=2.5.2
RUN conda clean -a -y

WORKDIR /data

CMD ["orthofinder"]
