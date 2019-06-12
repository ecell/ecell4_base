Installation
------------

### Windows

- Install Miniconda with **Python3.7 64-bit** from http://conda.pydata.org/miniconda.html
- Run the following commands in the command prompt.

    ```shell
    conda update conda python pip
    conda uninstall hdf5
    conda clean -a
    conda install hdf5 notebook
    pip install ecell4
    ```

If you need to create movie with E-Cell4, please install [ffmpeg windows build](http://ffmpeg.zeranoe.com/builds/) and add its path to your **USER** PATH enviromental variable.

### Mac

- Install Miniconda with **Python3.7 64-bit** from http://conda.pydata.org/miniconda.html
- Run the following commands in the Terminal app.

    ```shell
    ~/miniconda3/bin/conda update conda python pip
    ~/miniconda3/bin/conda install notebook
    ~/miniconda3/bin/pip install ecell4
    # If you need to create movie, install ffmpeg with homebrew
    brew install ffmpeg
    ```

### Any Linux (almost)

You do not need to install Miniconda on Linux.

```shell
python3 -m pip install ecell4
# If you need to create movie[, and use Ubuntu]
apt install ffmpeg
```
