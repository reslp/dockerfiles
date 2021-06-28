FROM reslp/mamba:0.7.0
MAINTAINER <philipp.resl@uni-graz.at>

# to fix locale warning in picard
RUN apt-get update && apt-get -y install locales
RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && \
    locale-gen
ENV LANG en_US.UTF-8  
ENV LANGUAGE en_US:en  
ENV LC_ALL en_US.UTF-8 


RUN mamba install picard=2.18.7


WORKDIR /data 

CMD ["picard"]
