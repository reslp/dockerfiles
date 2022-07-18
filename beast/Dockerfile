FROM ubuntu:20.04
MAINTAINER <philipp.resl@uni-graz.at>

ENV DEBIAN_FRONTEND=noninteractive
RUN apt update && apt install -y wget default-jre
WORKDIR /software
RUN wget https://github.com/CompEvol/beast2/releases/download/v2.6.7/BEAST.v2.6.7.Linux.tgz && tar xfzv BEAST.v2.6.7.Linux.tgz
ENV PATH=/software/beast/bin:$PATH

WORKDIR /data

CMD ["beast" , "-h"]
