FROM ubuntu:20.04

RUN apt-get update && apt-get install -y wget 

WORKDIR /software

RUN wget https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/table2asn_GFF/linux64.table2asn_GFF.gz
RUN gunzip linux64.table2asn_GFF.gz && chmod +x linux64.table2asn_GFF && mv linux64.table2asn_GFF /usr/local/bin

ENV PATH="/usr/local/bin/:$PATH"

WORKDIR /data

CMD ["linux64.table2asn_GFF", "-h"]
