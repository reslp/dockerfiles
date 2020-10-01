# Docker file for iqtree
# written by Philipp Resl, Dez. 2017
# last change: Oct. 2020
# with inspiration from: https://hub.docker.com/r/asmmhossain/pipelign/~/dockerfile/

# Use phusion/baseimage as base image
FROM phusion/baseimage:0.9.19

# Use baseimage-docker's init system.
CMD ["/sbin/my_init"]

# install necessary software
RUN apt-get update -qq && \
  	DEBIAN_FRONTEND=noninteractive apt-get install -yq --no-install-recommends \
 	build-essential \
 	python3-dev \
	git \
	wget \
	cmake && \
  apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Now Install IQ-TREE
WORKDIR /software
RUN wget https://github.com/Cibiv/IQ-TREE/releases/download/v2.0.7/iqtree-2.0.7-Linux.tar.gz && \
	tar xvfz iqtree-2.0.7-Linux.tar.gz && \
	cp iqtree-2.0.7-Linux/bin/iqtree2 /usr/local/bin && \
	ln -s /usr/local/bin/iqtree2 /usr/local/bin/iqtree
WORKDIR /data
# set entry point
ENTRYPOINT ["/usr/local/bin/iqtree"]
CMD ["-h"]
