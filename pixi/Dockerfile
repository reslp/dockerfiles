FROM ubuntu:26.04

RUN apt update -y && apt upgrade -y
RUN apt install -y curl

RUN curl -fsSL https://pixi.sh/install.sh | sh
ENV PATH="$PATH:/root/.pixi/bin"
