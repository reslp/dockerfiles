FROM ubuntu:20.04
MAINTAINER <philipp.resl@uni-graz.at>

# necessary for tzdata (installed as dependency of libopenmpi)
ENV TZ=Europe/Vienna
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y wget build-essential libopenmpi-dev

WORKDIR /software

RUN wget https://github.com/bayesiancook/pbmpi/archive/v1.8b.tar.gz && tar xvfz v1.8b.tar.gz
RUN cd pbmpi-1.8b/sources && make && cp -r ../data/* /usr/local/bin 

WORKDIR /data

CMD ["pb_mpi"]