FROM linuxbrew/brew:3.2.5
MAINTAINER <philipp.resl@uni-graz.at>

RUN brew update
RUN brew extract --version=1.6.6 nim homebrew/cask
RUN brew install nim@1.6.6

CMD ["nim"]
