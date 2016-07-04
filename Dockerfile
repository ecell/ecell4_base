FROM jupyter/minimal-notebook

RUN conda install -y matplotlib; pip install ecell; wget https://github.com/ecell/ecell4-notebooks/archive/master.zip; unzip master.zip; rm master.zip
