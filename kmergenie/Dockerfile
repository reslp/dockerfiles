FROM reslp/mamba:0.5.3
MAINTAINER <philipp.resl@uni-graz.at>

RUN mamba install -y python=2.7 r-base=3.5

WORKDIR /software
RUN apt-get update && apt-get install -y build-essential zlib1g-dev
RUN wget http://kmergenie.bx.psu.edu/kmergenie-1.7051.tar.gz && tar xvfz kmergenie-1.7051.tar.gz && cd kmergenie-1.7051 && make && rm ../kmergenie-1.7051.tar.gz

ENV PATH=/software/kmergenie-1.7051:$PATH
 
WORKDIR /data

CMD ["/bin/bash"]
