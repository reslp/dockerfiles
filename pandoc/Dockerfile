FROM reslp/mamba:0.21.1
MAINTAINER <philipp.resl@uni-graz.at>

RUN mamba install -c conda-forge pandoc=2.18
RUN DEBIAN_FRONTEND=noninteractive && apt update && apt install -y texlive
RUN apt install texlive-plain-generic
WORKDIR /data


CMD ["pandoc", "-h"]
