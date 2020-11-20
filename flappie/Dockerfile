FROM ubuntu:16.04


RUN apt-get update && apt-get install -y curl cmake libcunit1-dev libhdf5-dev libopenblas-dev
RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash && apt-get install -y git-lfs
RUN git clone https://github.com/nanoporetech/flappie
RUN cd flappie && git checkout 4de542f && make flappie

#add ont-fast4-api for convenience 
RUN apt-get install -y python3-pip
RUN pip3 install h5py==2.6 numpy==1.11 six==1.10.0 progressbar33==2.4 ont-fast5-api==3.1.6

ENV PATH="/flappie:$PATH"

WORKDIR /data
CMD ["/bin/bash"]
