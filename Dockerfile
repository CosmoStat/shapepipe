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
    apt-get install xterm -y && \
    apt-get install cmake protobuf-compiler -y && \
    apt-get clean

RUN apt-get install acl -y && \
    apt-get install sssd -y
ADD nsswitch.conf /etc/

RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && \
    locale-gen
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

FROM python:3.10-buster

RUN pip install poetry==1.7.0

ENV POETRY_NO_INTERACTION=1 \
    POETRY_VIRTUALENVS_IN_PROJECT=1 \
    POETRY_VIRTUALENVS_CREATE=1 \
    POETRY_CACHE_DIR=/tmp/poetry_cache

COPY . /

COPY pyproject.toml poetry.lock README.rst ./
RUN touch README.md

RUN poetry install --without dev --no-root && rm -rf $POETRY_CACHE_DIR

COPY shapepipe ./

RUN poetry install --without dev
