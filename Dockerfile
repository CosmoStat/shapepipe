FROM images.canfar.net/skaha/astroml:latest

# Metadata
LABEL maintainer="martin.kilbinger@cea.fr"
LABEL description="ShapePipe base image with common dependencies"

# Install system dependencies needed for ShapePipe and WeightWatcher
RUN echo "nameserver 8.8.8.8" > /etc/resolv.conf && \
    apt-get update -o Acquire::ForceIPv4=true -y --quiet && \
    apt-get install -y --no-install-recommends \
        psfex source-extractor \
        libproj-dev proj-bin && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Build and install WeightWatcher from source
ARG WW_VERSION=1.12
RUN cd /tmp && \
    wget --no-check-certificate https://github.com/astromatic/weightwatcher/archive/refs/tags/${WW_VERSION}.tar.gz && \
    tar -xzf ${WW_VERSION}.tar.gz && \
    rm ${WW_VERSION}.tar.gz
RUN cd /tmp/weightwatcher-${WW_VERSION} && \
    sed -i 's/^  prefstruct\tprefs;/extern prefstruct\tprefs;/' src/prefs.h && \
    sed -i 's/^char\t\tgstr\[MAXCHAR\];/extern char\t\tgstr[MAXCHAR];/' src/globals.h && \
    sed -i 's/^int\t\tbswapflag;/extern int\t\tbswapflag;/' src/fits/fitscat.h && \
    sed -i '/preflist\.h/a prefstruct\tprefs;' src/prefs.c && \
    sed -i '/xml\.h/a char\t\tgstr[MAXCHAR];' src/main.c && \
    sed -i '/fitscat\.h/a int\t\tbswapflag;' src/fits/fitscat.c && \
    ./configure --quiet && \
    make --quiet && \
    make install

# Install Python dependencies (including Jupyter Lab if desired)
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir \
        astropy==6.1.0 \
        cs_util==0.1.9 \
        galsim==2.5.3 \
        ipython==8.18.1 \
        joblib==1.4.2 \
        jupyterlab==4.3.1 \
        matplotlib==3.8.4 \
        mccd==1.2.4 \
        modopt==1.6.1 \
        mpi4py==4.0.3 \
        numpy==1.26.4 \
        numpydoc==1.2 \
        pandas==2.2 \
        pytest==8.3.3 \
        pytest-cov==5.0.0 \
        pytest-pycodestyle==2.4.1 \
        pytest-pydocstyle==2.4.0 \
        reproject==0.14.1 \
        sf_tools==2.0.4 \
        sip_tpv==1.1 \
        skaha==1.7.0 \
        sqlitedict==2.0.0 \
        termcolor==1.1.0 \
        tqdm==4.63.0 \
        treecorr==5.1.1 \
        vos==3.6.1.1 \
        "setuptools<81" \
        git+https://github.com/aguinot/ngmix@stable_version \
        git+https://github.com/tobias-liaudat/Stile@v0.1 \
        astroquery \
        canfar \
        h5py \
        PyQt5 \
        pyqtgraph \
        python-pysap \
        numba \
        fitsio

# Set working directory and copy source code
WORKDIR /app
COPY . /app/.
RUN chown -R root:root /app && chmod -R u+rwX /app

# Install ShapePipe in editable mode and expose scripts
RUN pip install --no-cache-dir -e . && \
    for ext in .py .sh .bash; do \
        for script in /app/scripts/*/*$ext; do \
            link_name=$(basename $script $ext); \
            ln -s $script /usr/local/bin/$link_name; \
        done; \
    done
