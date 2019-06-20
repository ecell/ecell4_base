FROM jupyter/minimal-notebook

USER root

RUN apt-get update && apt-get install -y ffmpeg unzip wget

USER jovyan

RUN conda install -y matplotlib && pip install ecell4 && \
    wget https://github.com/ecell/ecell4-base/archive/master.zip && \
    unzip master.zip
