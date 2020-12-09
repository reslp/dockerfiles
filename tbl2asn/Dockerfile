FROM ubuntu:20.04

RUN apt-get update && apt-get install -y wget libidn11

WORKDIR /software

RUN wget https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz
RUN gunzip linux64.tbl2asn.gz && chmod +x linux64.tbl2asn && mv linux64.tbl2asn /usr/local/bin/tbl2asn

ENV PATH="/usr/local/bin/:$PATH"

WORKDIR /data

CMD ["tbl2asn", "--help"]
