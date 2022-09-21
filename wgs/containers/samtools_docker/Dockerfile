FROM ubuntu:latest

RUN apt-get update && apt-get install --no-install-recommends -y \
 wget \
 bzip2 \
 zlib1g \
 zlib1g-dev \
 libncurses5-dev libncursesw5-dev \
 libbz2-dev \
 liblzma-dev \
 gcc\
 autoconf\
 make

RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 --no-check-certificate
RUN tar -xvf samtools-1.16.1.tar.bz2
RUN rm samtools-1.16.1.tar.bz2
WORKDIR samtools-1.16.1
RUN ./configure
RUN make
RUN make install

WORKDIR bwa
RUN  apt-get -y install bwa
