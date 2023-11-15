FROM continuumio/miniconda3

LABEL Description="ShapePipe Docker Image"
ENV SHELL /bin/bash

COPY ./environment.yml ./
RUN touch ./README.md

RUN conda env create --file environment.yml
