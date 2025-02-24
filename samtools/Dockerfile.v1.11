FROM reslp/mamba:0.5.3
MAINTAINER <philipp.resl@uni-graz.at>

ENV PATH /opt/conda/bin:$PATH

RUN mamba install -c bioconda -y samtools=1.11
	
WORKDIR /data

CMD ["samtools"]