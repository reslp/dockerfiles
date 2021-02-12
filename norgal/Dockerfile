FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
	apt-get install -y python3 python3-pip openjdk-8-jre && \
	ln -s /usr/bin/python3 /usr/bin/python
RUN python -m pip install -U matplotlib

WORKDIR /software

RUN apt-get install -y git && \
	git clone https://bitbucket.org/kosaidtu/norgal.git

RUN ln -s /software/norgal/norgal.py /usr/bin/norgal.py && \
	ln -s /software/norgal/norgal.py /usr/bin/norgal

WORKDIR /data
	
CMD ["/bin/bash"] 
