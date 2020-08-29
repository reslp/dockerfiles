FROM ubuntu:20.04
MAINTAINER <philipp.resl@uni-graz.at>

RUN apt-get update && apt-get install -y wget

WORKDIR /software

RUN wget https://megasun.bch.umontreal.ca/People/lartillot/www/phylobayes4.1c.tar.gz && tar xvfz phylobayes4.1c.tar.gz
RUN cp -r phylobayes4.1c/data/* /usr/bin && rm -r /software/*

WORKDIR /data

CMD ["pb"]