FROM ubuntu:latest

RUN apt-get update
RUN apt install -y default-jre
RUN apt-get install default-jdk
RUN apt-get -y install gzip
RUN apt-get install zip
RUN apt-get install unzip
RUN apt-get install wget
WORKDIR snpEff
RUN wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
RUN unzip snpEff_latest_core.zip
RUN rm snpEff_latest_core.zip
