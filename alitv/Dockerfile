FROM ubuntu:20.04
MAINTAINER <philipp.resl@uni-graz.at>

ENV DEBIAN_FRONTEND=noninteractive
RUN apt update -y
RUN apt install -y libyaml-perl libhash-merge-perl bioperl git
RUN apt install -y cpanminus build-essential
WORKDIR /software
RUN git clone --recursive https://github.com/AliTVTeam/AliTV
RUN cd AliTV && git checkout d06bd45 && cd /software/AliTV/AliTV-perl-interface/ && cpanm --installdeps . --force
ENV PATH="$PATH:/software/AliTV/AliTV-perl-interface/bin"
WORKDIR /data

CMD ["alitv.pl"]
