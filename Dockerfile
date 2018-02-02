FROM jupyter/minimal-notebook

USER root

RUN apt-get update && apt-get install -y ffmpeg

USER jovyan

RUN conda install -y matplotlib && pip install ecell && \
    mkdir notebooks && cd notebooks && git init && \
    git remote add origin -f https://github.com/ecell/ecell4.git && \
    git config core.sparsecheckout true && \
    echo "readthedocs/examples/*" >> .git/info/sparse-checkout && \
    echo "readthedocs/tutorials/*" >> .git/info/sparse-checkout && \
    git pull origin master
