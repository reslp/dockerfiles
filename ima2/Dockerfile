FROM ubuntu:20.04
MAINTAINER <philipp.resl@uni-graz.at>

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y wget build-essential

WORKDIR /software
RUN wget https://github.com/jodyhey/archived/raw/master/IMA2/ima2-8.27.12.tar.gz && tar xvfz ima2-8.27.12.tar.gz 
RUN cd ima2-8.27.12 && ./configure && make
ENV PATH=$PATH:"/software/ima2-8.27.12/src"

WORKDIR /data

CMD ["IMa2"]
