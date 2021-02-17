FROM reslp/mamba:0.7.0
MAINTAINER <philipp.res@uni-graz.at>

RUN mamba install -y -c bioconda novoplasty=4.2
WORKDIR /data

CMD ["NOVOPlasty.pl"]
