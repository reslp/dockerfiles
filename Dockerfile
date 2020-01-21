FROM ubuntu:18.04

MAINTAINER <resl.philipp@bio.lmu.de>

WORKDIR /data

#these two layers should take care of all the python and perl dependencies:
RUN apt-get update && apt-get -y upgrade && apt-get -y install wget && apt-get -y install python2.7 python-pip && pip install funannotate
RUN apt-get install -y bioperl && apt-get install -y cpanminus && cpanm File::Which Hash::Merge JSON Logger::Simple Parallel::ForkManager Scalar::Util::Numeric Text::Soundex

#Software dependencies:
#CodingQuarry
RUN wget https://sourceforge.net/projects/codingquarry/files/CodingQuarry_v2.0.tar.gz && tar xvfz CodingQuarry_v2.0.tar.gz && cd CodingQuarry_v2.0 && make 
ENV QUARRY_PATH="/data/CodingQuarry_v2.0/QuarryFiles"
ENV PATH="/data/CodingQuarry_v2.0:$PATH"

#Trinity
RUN apt-get -y install cmake && apt-get -y install zlib1g-dev && wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/v2.8.6/trinityrnaseq-v2.8.6.FULL.tar.gz && tar xvfz trinityrnaseq-v2.8.6.FULL.tar.gz && cd trinityrnaseq-v2.8.6 && make
ENV TRINITY_HOME="/data/trinityrnaseq-v2.8.6"
ENV PATH="/data/trinityrnaseq-v2.8.6:$PATH"

ENTRYPOINT ["/bin/bash"]