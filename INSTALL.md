Installation
------------

### Requirements

#### Minimum requirements

**Note that we do NOT support Python2.7 and 32bit versions.**

- Python (3.5, 3.6, 3.7 are supported [3.5 is supported only on Linux])
- pip (8.1 or later)
- HDF5 (required only on **Windows**.)

#### Optional requirements

We strongly recommend that you run E-Cell4 with [Jupyter Notebook](http://jupyter.org/).
Some E-Cell4 functions (for movie, datastore) optionally depend on

- ffmpeg
- pandas

### Windows

We recommend that you install [Miniconda](http://conda.pydata.org/miniconda.html) to manage Python packages.

- Install Miniconda for Windows from http://conda.pydata.org/miniconda.html
- Run the following commands on command prompt

    ```shell
    conda uninstall hdf5
    conda clean -a
    conda install hdf5 notebook
    pip install ecell4
    ```

If you need to create movie with E-Cell4, please install [ffmpeg windows build](http://ffmpeg.zeranoe.com/builds/) and add its path to your **USER** PATH enviromental variable.

### Mac

We recommend that you install [Miniconda](http://conda.pydata.org/miniconda.html) to manage Python packages.
After installing Miniconda, run the following commands in your terminal.

```shell
# After installing Miniconda2 or Miniconda3 (Here we assume that you installed Miniconda3).
~/miniconda3/bin/conda install notebook
~/miniconda3/bin/pip install ecell4
# If you need to create movie, install ffmpeg with homebrew
brew install ffmpeg
```

### Linux

You do not need to use conda in Linux.

```shell
python3 -m pip install ecell4
# If you use Ubuntu
apt install ffmpeg
```
