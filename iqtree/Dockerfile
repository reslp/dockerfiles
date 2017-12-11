# Use phusion/baseimage as base image
FROM phusion/baseimage:0.9.19

	
# Use baseimage-docker's init system.
CMD ["/sbin/my_init"]


RUN apt-get update -qq && \
  	DEBIAN_FRONTEND=noninteractive apt-get install -yq --no-install-recommends \
 	build-essential \
	#default-jdk \
	#perl \
 	python3-dev \
	git \
	wget \
	cmake
	
RUN mkdir /home/docker
RUN mkdir /home/docker/programs

# Install IQ-TREE
RUN git clone https://github.com/Cibiv/IQ-TREE.git && \
	cd IQ-TREE && \
	mkdir build && \
	cd build && \
	cmake -DIQTREE_FLAGS=omp .. && \
	make && \
	cp iqtree-omp /usr/local/bin && \
	rm -rf IQ-TREE

RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENTRYPOINT ["/usr/local/bin/iqtree-omp"]
CMD ["-h"]