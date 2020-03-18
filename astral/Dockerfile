FROM ubuntu:18.04

RUN apt-get update && \
	apt-get install -y --no-install-recommends wget unzip default-jre && \
	wget https://github.com/smirarab/ASTRAL/archive/v5.7.1.tar.gz && \
	tar xvfz v5.7.1.tar.gz && \
	cd ASTRAL-5.7.1 && \
	unzip Astral.5.7.1.zip

CMD ["java", "-jar", "/ASTRAL-5.7.1/Astral/astral.5.7.1.jar"]