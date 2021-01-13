FROM ubuntu:20.04
MAINTAINER <philipp.resl@uni-graz.at>

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y build-essential default-jre python3 && \
	apt-get clean && apt-get purge && \
        rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
