FROM ubuntu:22.04

RUN apt update && apt install -y wget

WORKDIR /software

RUN apt install -y python3 python3.10-venv build-essential python3-dev
RUN wget -O oma.tgz https://omabrowser.org/standalone/OMA.2.5.0.tgz
RUN tar xvzf oma.tgz
RUN cd OMA.2.5.0 && ./install.sh # will install to /usr/local/OMA
ENV PATH="$PATH:/usr/local/OMA/bin"

WORKDIR /data

CMD ["OMA"]


