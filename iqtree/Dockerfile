# Docker file for IQ-tree
# written by Philipp Resl, Dez. 2017
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
	cmake
	
# Now Install IQ-TREE
RUN git clone https://github.com/Cibiv/IQ-TREE.git && \
	cd IQ-TREE && \
	mkdir build && \
	cd build && \
	cmake -DIQTREE_FLAGS=omp .. && \
	make && \
	cp iqtree-omp /usr/local/bin && \
	rm -rf IQ-TREE

# clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# set entry point
ENTRYPOINT ["/usr/local/bin/iqtree-omp"]
CMD ["-h"]