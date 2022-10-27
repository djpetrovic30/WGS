FROM ubuntu:latest

RUN export DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt install -y default-jre
RUN apt-get install -y default-jdk
RUN apt-get -y install gzip
RUN apt-get -y install zip
RUN apt-get -y install unzip
RUN apt-get -y install wget
WORKDIR ./opt
RUN wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
RUN unzip snpEff_latest_core.zip
RUN rm snpEff_latest_core.zip
WORKDIR ./snpEff
RUN java -jar snpEff.jar download -v GRCh37.75
