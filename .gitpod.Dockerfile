# Image source code: https://github.com/axonasif/workspace-images/tree/tmp
# Also see https://github.com/gitpod-io/workspace-images/issues/1071
FROM axonasif/workspace-base@sha256:8c057b1d13bdfe8c279c68aef8242d32110c8d5310f9a393f9c0417bc61367d9

# Set user
USER gitpod

# Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py311_24.1.2-0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p $HOME/miniconda && \
    rm ~/miniconda.sh

# Add Miniconda to PATH
ENV PATH="$HOME/miniconda/bin:$PATH"

# Initialize conda in bash config files:
RUN conda init bash

# Set up Conda channels
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --set channel_priority strict

# Set libmamba as solver
RUN conda config --set solver libmamba

# Persist ~/ (HOME) and lib
RUN echo 'create-overlay $HOME /lib' > "$HOME/.runonce/1-home-lib_persist"

# Referenced in `.vscode/settings.json`
ENV PYTHON_INTERPRETER="$HOME/miniconda/bin/python"
# Pycharm recognizes this variables
ENV PYCHARM_PYTHON_PATH="${PYTHON_INTERPRETER}"