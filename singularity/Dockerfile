FROM centos:7
MAINTAINER <resl.philipp@bio.lmu.de>

RUN yum groupinstall -y 'Development Tools' && \
	yum install -y epel-release && \
	yum install -y golang openssl-devel libuuid-devel libseccomp-devel squashfs-tools cryptsetup

RUN mkdir -p ${GOPATH}/src/github.com/sylabs && \
	cd ${GOPATH}/src/github.com/sylabs && \
	git clone https://github.com/sylabs/singularity.git && \
	cd singularity && \
	git checkout v3.6.3 && \
	./mconfig && \
    make -C ./builddir && \
    make -C ./builddir install

ENV HOME="/data"
WORKDIR /data

ENTRYPOINT ["singularity"]