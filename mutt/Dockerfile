FROM ubuntu:20.04

RUN apt update && apt install -y build-essential wget libsqlite3-dev libgpgme-dev libncurses5-dev libncursesw5-dev libssl-dev ca-certificates gnutls-bin
RUN apt install -y libdb5.3-dev libsasl2-dev 
RUN mkdir /software && cd /software && wget http://ftp.mutt.org/pub/mutt/mutt-2.2.6.tar.gz && tar xvfz mutt-2.2.6.tar.gz && cd mutt-2.2.6
RUN cd /software/mutt-2.2.6 && ./configure --with-sqlite3 --enable-autocrypt --enable-gpgme --enable-pop --enable-imap --enable-smtp --enable-compressed --enable-sidebar --with-ssl=/usr/bin/ --enable-hcache --with-sasl && make install

# set locale:
RUN apt install locales locales-all
ENV LC_ALL de_AT.UTF-8
ENV LANG de_AT.UTF-8
ENV LANGUAGE de_AT.UTF-8

RUN apt install -y vim gnupg
