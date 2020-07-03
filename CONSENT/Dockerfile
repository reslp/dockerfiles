FROM ubuntu:18.04

RUN apt-get update && apt-get install -y git build-essential zlib1g-dev cmake
RUN git clone --recursive https://github.com/morispi/CONSENT
RUN cd /CONSENT/ && git checkout v2.1 && ./install.sh
ENV PATH="/CONSENT:$PATH"
ENV PATH="/CONSENT/bin:$PATH"
ENV PATH="/CONSENT/minimap2:$PATH"
ENV PATH="/CONSENT/BMEAN:$PATH"

CMD ["CONSENT-correct"]
