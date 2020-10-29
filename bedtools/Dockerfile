FROM ubuntu:20.04

RUN apt-get update && apt-get install -y wget

WORKDIR /software

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools.static.binary && cp bedtools.static.binary /usr/bin/bedtools && chmod +x /usr/bin/bedtools

WORKDIR /data

CMD ["bedtools"]
