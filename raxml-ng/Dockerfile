FROM ubuntu:18.04

RUN apt-get update && \
	apt-get install -y flex bison libgmp3-dev build-essential git wget unzip
WORKDIR /software

RUN wget https://github.com/amkozlov/raxml-ng/releases/download/1.1.0/raxml-ng_v1.1.0_linux_x86_64.zip
RUN unzip raxml-ng_v1.1.0_linux_x86_64.zip

ENV PATH="/software:$PATH"

WORKDIR /data
CMD ["raxml-ng"]
