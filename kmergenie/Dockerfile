FROM reslp/mamba:0.5.3
MAINTAINER <philipp.resl@uni-graz.at>

RUN mamba install -y -c bioconda kmergenie=1.7051 && \
	mamba clean -a -y

WORKDIR /data

CMD ["kmergenie"]
