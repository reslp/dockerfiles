FROM ubuntu:18.04
# this follows what is described in an issue on the bonito github page, with some minor adjustments:
# https://github.com/nanoporetech/bonito/issues/70

USER root

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      python3-minimal \
      python3-pip \
      git \
      build-essential \
      python3-dev \
      libz-dev

RUN git clone -b openvino --depth 1 https://github.com/dkurt/bonito /software/bonito

RUN cd /software/bonito && \
    python3 -m pip install --upgrade pip && \
    python3 -m pip install setuptools && \
    python3 -m pip install -r requirements.txt && \
    python3 setup.py develop && \
    bonito download --models
	
RUN apt-get install -y wget libssl-dev

WORKDIR /software

RUN wget https://github.com/Kitware/CMake/releases/download/v3.18.5/cmake-3.18.5.tar.gz && \
	tar xf cmake-3.18.5.tar.gz && \
	cd cmake-3.18.5 && \
	./configure && \
	make && \
	make install

#recompile openvino to be able to use v3 models also for basecalling:
RUN cd /software && git clone --depth 1 https://github.com/openvinotoolkit/openvino/ && \
	cd openvino && \
	git submodule update --init --recursive && \
	mkdir build && cd build && \
	cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_CLDNN=OFF -DENABLE_VPU=OFF -DENABLE_TESTS=OFF -DENABLE_PYTHON=ON -DENABLE_OPENCV=OFF -DENABLE_SAMPLES=OFF -DPYTHON_EXECUTABLE=`which python3` -DPYTHON_LIBRARY=/usr/lib/x86_64-linux-gnu/libpython3.6m.so -DPYTHON_INCLUDE_DIR=/usr/include/python3.6 .. && \
	make -j$(nproc)
	


RUN useradd -ms /bin/bash -G users bonito && \
    chown bonito -R /home/bonito

#USER bonito

ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
ENV LD_LIBRARY_PATH=/software/openvino/bin/intel64/Release/lib/:$LD_LIBRARY_PATH
ENV PYTHONPATH=/software/openvino/bin/intel64/Release/lib/python_api/python3.6/:$PYTHONPATH


WORKDIR /data

CMD ["bonito"]

