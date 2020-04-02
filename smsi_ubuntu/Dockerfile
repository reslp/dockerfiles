#FROM ubuntu:18.04
FROM continuumio/miniconda3:4.7.10
MAINTAINER <philipp.resl@uni-graz.at>

RUN apt-get update && apt-get install -y \
    build-essential \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    squashfs-tools \
    libseccomp-dev \
    pkg-config \
	wget \
	git \
	cryptsetup-bin

RUN wget https://dl.google.com/go/go1.13.1.linux-amd64.tar.gz && \
	tar -xvf go1.13.1.linux-amd64.tar.gz && \
	mv go /usr/local
	
ENV GOROOT="/usr/local/go"
ENV GOPATH=$HOME/go
ENV PATH="$GOPATH/bin:$GOROOT/bin:$PATH"


RUN mkdir -p ${GOPATH}/src/github.com/sylabs && \
	cd ${GOPATH}/src/github.com/sylabs && \
	git clone https://github.com/sylabs/singularity.git && \
	cd singularity && \
	git checkout v3.5.2 && \
	./mconfig && \
    make -C ./builddir && \
    make -C ./builddir install

ENV HOME="/data"
WORKDIR /data

RUN conda config --set notify_outdated_conda false
RUN conda install -c bioconda -c conda-forge snakemake=5.13.0

CMD ["/bin/bash"]