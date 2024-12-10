FROM python:3.9-bookworm

LABEL Description="Conda-Free ShapePipe Docker Image"
ENV SHELL=/bin/bash

# Install system dependencies
RUN apt-get update -y --quiet --fix-missing && \
    apt-get dist-upgrade -y --quiet --fix-missing && \
    apt-cache policy autconf && \ 
    apt-get install -y --quiet \
    apt-utils \
    autoconf \
    automake \
    build-essential \
    cmake \
    curl \
    ffmpeg \
    g++ \
    gcc  \
    gfortran \
    git-lfs \
    libatlas-base-dev \
    libblas-dev \
    liblapack-dev \
    libcfitsio-dev \
    libfftw3-bin \
    libfftw3-dev \
    libgl1-mesa-glx \
    libtool \
    libtool-bin \
    libtool-doc \
    locales \
    locate \
    make \
    openmpi-bin \
    libopenmpi-dev \
    pkg-config \
    protobuf-compiler \
    psfex=3.21.1-1 \
    source-extractor=2.25.0+ds-3 \
    weightwatcher=1.12+dfsg-3 \
    vim \
    xterm && \
    apt-get clean -y && \
    apt-get autoremove --purge --quiet -y && \
    rm -rf /var/lib/apt/lists/* /var/tmp/*

# Install CDS client by hand
RUN cd /tmp && \
    curl -O http://cdsarc.u-strasbg.fr/ftp/pub/sw/cdsclient.tar.gz && \
    tar xvfz cdsclient.tar.gz \ 
    && cd cdsclient-* \
    && ./configure && make && make install \
    && rm -rf /tmp/*

# Install python dependencies
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir \
    astropy==5.2 \
    cs_util==0.1.0 \
    galsim==2.2.6 \
    ipython==8.18.1 \
    joblib==1.1.0 \
    jupyterlab==4.3.1 \
    matplotlib==3.8.4 \
    mccd==1.2.4 \
    modopt==1.6.1 \
    mpi4py==3.1.3 \
    numba==0.58.1 \
    numpy==1.26.4 \
    numpydoc==1.2  \
    pandas==1.4.1  \
    PyQt5==5.15.6 \
    pyqtgraph==0.12.4 \
    pytest==8.3.3 \
    pytest-cov==5.0.0 \
    pytest-pycodestyle==2.4.1 \
    pytest-pydocstyle==2.4.0 \
    python-pysap==0.2.1 \
    reproject==0.8 \
    sf_tools==2.0.4 \
    sip_tpv==1.1 \
    skaha==1.4.3 \
    sqlitedict==2.0.0 \
    termcolor==1.1.0 \
    tqdm==4.63.0 \
    treecorr==4.2.6 \
    vos==3.6.1.1 \
    git+https://github.com/aguinot/ngmix@stable_version \
    git+https://github.com/tobias-liaudat/Stile@v0.1


WORKDIR /app
COPY . /app/.

# Install shapepipe and symlink scripts
RUN pip install --no-cache-dir -e . && \ 
    for ext in .py .sh .bash; do \
    for script in /app/scripts/*/*$ext; do \
    link_name=`basename $script $ext`; \
    ln -s $script /usr/local/bin/$link_name; \
    done; \
    done