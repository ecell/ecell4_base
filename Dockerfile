FROM cytoscape/jupyter-cytoscape

USER root

RUN apt-get update && apt-get install -y libav-tools

USER $NB_USER

RUN conda install -y matplotlib && pip install ecell
