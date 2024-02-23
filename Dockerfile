FROM continuumio/miniconda3

LABEL Description="ShapePipe Docker Image"
ENV SHELL /bin/bash

ARG CC=gcc-9
ARG CXX=g++-9

# gcc < 10 is required to compile ww
ENV CC=gcc-9
ENV CXX=g++-9

RUN apt-get update --allow-releaseinfo-change && \
    apt-get update && \
    apt-get upgrade -y && \
    apt-get install apt-utils -y && \
    apt-get install make -y && \
    apt-get install automake -y && \
    apt-get install autoconf -y && \
    apt-get install gcc-9 g++-9 -y && \
    apt-get install locales -y && \
    apt-get install libgl1-mesa-glx -y && \
    apt-get install xterm -y && \
    apt-get install cmake protobuf-compiler -y && \
    apt-get install libtool libtool-bin libtool-doc -y && \
    apt-get install libfftw3-bin libfftw3-dev -y && \
    apt-get install libatlas-base-dev liblapack-dev libblas-dev -y && \
    apt-get install vim -y && \
    apt-get install locate -y && \
    apt-get install curl -y && \
    apt-get install acl -y && \
    apt-get install sssd -y && \
    apt-get clean

ADD nsswitch.conf /etc/

RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && \
    locale-gen
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

SHELL ["/bin/bash", "--login", "-c"]

COPY ./environment.yml ./
COPY install_shapepipe README.rst setup.py setup.cfg ./
RUN touch ./README.md

RUN conda update -n base -c defaults conda -c defaults
RUN conda env create --file environment.yml

COPY shapepipe ./shapepipe
COPY scripts ./scripts

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "shapepipe", "/bin/bash", "-c"]
RUN pip install jupyter
