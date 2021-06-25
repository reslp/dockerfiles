FROM ubuntu:18.04

RUN apt-get update -y && \
	apt-get upgrade -y && \
	apt-get install -y python-pip && \
	pip install speedtest-cli==2.1.3
	
CMD ["speedtest-cli"]
