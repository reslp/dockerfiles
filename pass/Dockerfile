FROM ubuntu:22.04

RUN apt update && apt install -y build-essential wget

RUN mkdir /software && cd /software && wget https://git.zx2c4.com/password-store/snapshot/password-store-1.7.4.tar.xz && tar xvf password-store-1.7.4.tar.xz && cd password-store-1.7.4 && make install

RUN apt install -y vim gnupg-agent pinentry-curses
RUN apt install -y tree
RUN apt clean

ENV HOME=/data

CMD ["pass"]
