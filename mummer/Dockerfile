FROM reslp/mamba:0.5.3
MAINTAINER <philipp.resl@uni-graz.at>

RUN mamba install -y -c bioconda mummer4=4.0.0beta2
RUN apt-get update && apt-get install -y gnuplot gnuplot-x11

# setup for a fake x server:
RUN apt-get install -y xvfb
ENV DISPLAY :99
ADD run.sh /run.sh 
RUN chmod a+x /run.sh

RUN apt-get install -y fig2dev xfig 

WORKDIR /data

CMD ["/bin/bash", "/run.sh"]
