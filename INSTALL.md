Installation
------------

### Requirements

#### Minimum requirements

- Python (2.7 and 3.4, 3.5, 3.6 both major versions are supported [3.4 is supported only on Linux, Mac does not support 3.4 and 3.5])
- pip (8.1 or later)
- HDF5 (1.8.17, required only on **Windows**.)

#### Optional requirements

We strongly recommend that you run E-Cell4 with [Jupyter Notebook](http://jupyter.org/).
Some E-Cell4 functions (for visualization, datastore) optionally depend on

- matplotlib (**1.5.1** or later)
- ffmpeg
- pandas

### Windows

We recommend that you install [Miniconda](http://conda.pydata.org/miniconda.html) to manage Python packages.
**Note that we do not support Python2.7 64bit for Windows.**

(**Important**) E-Cell4 for Windows needs `HDF5` version **1.8.17**.

- Install Miniconda for Windows from http://conda.pydata.org/miniconda.html
- Run the following commands on **Anaconda Prompt** (You can run Anaconda Prompt from Windows Start Menu).

    ```shell
    conda uninstall hdf5
    conda clean -a
    conda install hdf5=1.8.17
    conda install matplotlib notebook
    pip install ecell
    ```

If you use animated visualization with E-Cell4, please install [ffmpeg windows build](http://ffmpeg.zeranoe.com/builds/) and add its path to your **USER** PATH enviromental variable.

### Mac or Linux

We recommend that you install [Miniconda](http://conda.pydata.org/miniconda.html) to manage Python packages.
After installing Miniconda, run the following commands in your terminal.

(NOTICE for Mac users) We do not provide **Python3.5 whl for Mac**. Instead we provide **Python3.6 whl for Mac**. To use Python3.6 enviroment, please refer to http://conda.pydata.org/docs/py2or3.html . Continuum.io already offers Python3.6 conda packages.

```shell
# After installing Miniconda2 or Miniconda3 (Here we assume that you installed Miniconda3).
~/miniconda3/bin/conda install matplotlib notebook

# Download E-Cell4 whl file for your Python version from https://github.com/ecell/ecell4/releases before running this command.
~/miniconda3/bin/pip install ecell

# If you want animation support, install ffmpeg with homebrew
brew install ffmpeg
# or if you use Ubuntu Linux
# apt install ffmpeg
```

### Docker

If you're familiar with Docker, the following commands should work in most cases:

```shell
docker pull ecell/ecell4
docker run -d -p 8888:8888 ecell/ecell4 start-notebook.sh --NotebookApp.token=''
```

and open a web browser to `http://localhost:8888` .

Our Docker image is based on **Minimal Jupyter Notebook Stack**. See https://github.com/jupyter/docker-stacks/tree/master/base-notebook or [Our Wiki page](https://github.com/ecell/ecell4/wiki/Security-in-the-Docker-Jupyter-notebook-server) for more details on the Docker command options.
