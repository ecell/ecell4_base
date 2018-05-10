FROM jupyter/minimal-notebook

USER root

RUN apt-get update && apt-get install -y ffmpeg

USER jovyan

RUN conda install -y matplotlib && pip install ecell
