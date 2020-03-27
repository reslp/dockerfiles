FROM ubuntu:18.04

RUN apt-get update -y && \
	apt-get upgrade -y && \
	apt-get install -y wget && \
	wget http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz && \
	tar xvfz pal2nal.v14.tar.gz && \
	mv /pal2nal.v14/pal2nal.pl /usr/local/bin
	
CMD ["pal2nal.pl"]
