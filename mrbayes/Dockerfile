# this creates a container with mrbayes-mpi installed
# can be called with:
# docker run -v <path_to_local_folder>:/input/ -it mb_ubuntu /sbin/my_init -- mpirun --allow-run-as-root -np 8 mb /input/<nexus_file.nex>

FROM phusion/baseimage:0.9.21

# Use baseimage-docker's init system.
CMD ["/sbin/my_init"]

RUN  apt-get update \
  && apt-get install -y wget \
  && apt-get install -y openmpi-bin \
  && apt-get install -y openmpi-doc \
  && apt-get install -y libopenmpi-dev \
  && apt-get install -y autoconf \
  && apt-get install -y build-essential \
  && wget https://sourceforge.net/projects/mrbayes/files/mrbayes/3.2.6/mrbayes-3.2.6.tar.gz \
  && tar -xf mrbayes-3.2.6.tar.gz \
  && cd mrbayes-3.2.6/src/ \
  && autoconf \
  && ./configure --enable-mpi=yes --with-beagle=no \
  && make \
  && make install

RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
