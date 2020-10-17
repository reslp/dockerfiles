FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt update && \
	apt install -y build-essential wget vim  
RUN apt install -y gdb

WORKDIR /data
CMD ["/bin/bash"]

