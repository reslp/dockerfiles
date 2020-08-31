FROM ubuntu:20.04
MAINTAINER <philipp.resl@uni-graz.at>

RUN apt-get update

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get install -y wget build-essential ruby-full cpanminus libexpat1-dev libssl-dev bioperl
#RUN apt-get install -y --no-install-recommends cpanminus 

WORKDIR /software

RUN wget http://guidance.tau.ac.il/guidance.v2.02.tar.gz && \
	tar xvfz guidance.v2.02.tar.gz && \
	cd guidance.v2.02 && \
	make
RUN wget https://mafft.cbrc.jp/alignment/software/mafft_7.471-1_amd64.deb && \
	dpkg -i mafft_7.471-1_amd64.deb

RUN cpanm Bio::Perl

WORKDIR /data

CMD ["perl", "/software/guidance.v2.02/www/Guidance/guidance.pl"]