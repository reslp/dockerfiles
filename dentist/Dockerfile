FROM reslp/mamba:0.15.3
MAINTAINER <philipp.resl@uni-graz.at>

RUN mamba install snakemake=6.0.2

RUN wget https://github.com/a-ludi/dentist/releases/download/v3.0.0/dentist.v3.0.0.x86_64.tar.gz

RUN tar -xzf dentist.v3.0.0.x86_64.tar.gz

ENV PATH="/dentist.v3.0.0.x86_64/bin:$PATH"

CMD ["dentist"]
