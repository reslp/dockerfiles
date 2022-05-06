FROM reslp/mamba:0.21.1

RUN mamba install python=3.7.13
RUN mamba install -c bioconda dendropy=4.5.2
RUN mamba install -c etetoolkit ete3=3.1.2
RUN mamba install pandas

WORKDIR /data

CMD ["python"]
