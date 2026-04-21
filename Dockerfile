FROM images.canfar.net/skaha/astroml:latest

# Metadata
LABEL maintainer="martin.kilbinger@cea.fr"
LABEL description="ShapePipe base image with common dependencies"

# Install system dependencies needed for ShapePipe and WeightWatcher
RUN apt-get update -o Acquire::ForceIPv4=true -y --quiet && \
    apt-get install -y --no-install-recommends \
    psfex source-extractor \
    libproj-dev proj-bin && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

#RUN echo "nameserver 8.8.8.8" > /etc/resolv.conf && \
    #apt-get update -o Acquire::ForceIPv4=true -y --quiet && \
    #apt-get install -y --no-install-recommends \
        #psfex source-extractor \
        #libproj-dev proj-bin && \
    #apt-get clean && rm -rf /var/lib/apt/lists/*

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

# Ensure astroml:latest conda Python 3.12 is used (Docker RUN does not source conda init)
ENV PATH /opt/conda/bin:$PATH

# Upgrade pip and install tools not part of the ShapePipe package
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir \
        ipython==8.18.1 \
        jupyterlab==4.3.1 \
        snakemake==8.27.1

# Set working directory and copy source code
WORKDIR /app
COPY . /app/.
RUN chown -R root:root /app && chmod -R u+rwX /app

# Install ShapePipe and its dependencies (including fitsio optional extra)
RUN pip install --no-cache-dir -e ".[fitsio]" && \
    for ext in .py .sh .bash; do \
        for script in /app/scripts/*/*$ext; do \
            link_name=$(basename $script $ext); \
            ln -s $script /usr/local/bin/$link_name; \
        done; \
    done
