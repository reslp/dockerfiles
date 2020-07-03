FROM ubuntu:18.04

RUN apt-get update && apt-get install -y git build-essential zlib1g-dev cmake
RUN git clone --recursive https://github.com/holtjma/fmlrc
RUN cd /fmlrc/ && git checkout v1.0.0 && make
ENV PATH="/fmlrc:$PATH"

CMD ["fmlrc"]
