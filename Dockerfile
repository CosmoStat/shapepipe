FROM continuumio/miniconda3

LABEL Description="ShapePipe Docker Image"
WORKDIR /home
ENV SHELL /bin/bash

ARG CC=gcc-9
ARG CXX=g++-9

RUN apt-get update --allow-releaseinfo-change && \
    apt-get update && \
    apt-get upgrade -y && \
    apt-get install make -y && \
    apt-get install gcc-9 g++-9 -y && \
    apt-get install locales -y && \
    apt install libgl1-mesa-glx -y && \
    apt-get clean

RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && \
    locale-gen
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

COPY . /home

RUN ./install_shapepipe --develop --vos
