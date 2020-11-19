FROM ubuntu:16.04


RUN apt-get update && apt-get install -y curl cmake libcunit1-dev libhdf5-dev libopenblas-dev
RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash && apt-get install -y git-lfs
RUN git clone https://github.com/nanoporetech/flappie
RUN cd flappie && git checkout 4de542f && make flappie

ENV PATH="/flappie:$PATH"

WORKDIR /data
CMD ["/bin/bash"]
