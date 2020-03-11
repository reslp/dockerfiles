FROM ubuntu:18.04
MAINTAINER <philipp.resl@uni-graz.at>

RUN apt-get update && \
	apt-get install -y wget gfortran build-essential --no-install-recommends && \
	wget http://ceiba.biosci.arizona.edu/r8s/r8s1.81.tar.gz && \
	tar xvfz r8s1.81.tar.gz && \
	cd /r8s1.81/src && \
	mkdir /usr/include/sys && \
	cp /usr/include/errno.h /usr/include/sys && \
	sed -i 's/continuousML.o:continuousML.h/#continuousML.o:continuousML.h/' Makefile.linux  && \
	sed -i 's/blas.o tn.o TNwrapper.o continuousML.o ancestral.o covarion.o/blas.o tn.o TNwrapper.o ancestral.o covarion.o #continuousML.o/' Makefile.linux  && \
	make -f Makefile.linux && \
	mv r8s /bin
	
CMD ["r8s"]