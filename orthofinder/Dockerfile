FROM reslp/mamba:1.4.1

MAINTAINER <philipp.resl@uni-graz.at>

RUN mamba install -c bioconda orthofinder=2.5.4
RUN conda clean -a -y

WORKDIR /data

CMD ["orthofinder"]
