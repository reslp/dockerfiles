FROM ubuntu:20.04

RUN apt update && apt install -y openssl ca-certificates xclip wget
RUN wget https://github.com/timvisee/ffsend/releases/download/v0.2.76/ffsend-v0.2.76-linux-x64-static && cp ffsend-v0.2.76-linux-x64-static /usr/bin/ffsend && chmod +x /usr/bin/ffsend
WORKDIR /data

CMD ["ffsend"]
