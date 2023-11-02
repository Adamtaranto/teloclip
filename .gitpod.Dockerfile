# Use the Gitpod base image for Python 3.10
FROM gitpod/workspace-python-3.12

# Set user
USER gitpod

# Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p $HOME/miniconda && \
    rm ~/miniconda.sh

# Add Miniconda to PATH
ENV PATH="$HOME/miniconda/bin:$PATH"

ENV PYTHONPATH="$HOME/miniconda/bin"

# Initialize conda in bash config fiiles:
RUN conda init bash

# Set up Conda channels
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --set channel_priority strict

# Set libmamba as solver
RUN conda config --set solver libmamba

# Remove the undesired location from PATH
RUN export PATH=$(echo $PATH | tr ':' '\n' | grep -v '/home/gitpod/.pyenv/shims' | tr '\n' ':')
