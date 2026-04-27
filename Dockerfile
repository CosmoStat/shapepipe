FROM python:3.12-slim-bookworm

# Metadata
LABEL maintainer="martin.kilbinger@cea.fr"
LABEL description="ShapePipe base image — slim Python + uv-frozen deps"

ENV SHELL=/bin/bash \
    QT_QPA_PLATFORM=offscreen \
    PIP_NO_CACHE_DIR=1 \
    DEBIAN_FRONTEND=noninteractive

# System dependencies. Three categories:
#  - astromatic binaries (psfex, source-extractor, weightwatcher) ship as
#    Debian packages on bookworm; preferred over building from source.
#  - compilers and dev libs needed to build the heavier wheels (galsim,
#    mpi4py, python-pysap, fitsio).
#  - libgl1, proj, fftw at runtime for skyproj/PyQt5/galsim.
RUN apt-get update -y --quiet && \
    apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        gfortran \
        git \
        wget \
        pkg-config \
        libfftw3-dev libfftw3-bin \
        libgsl-dev \
        libcfitsio-dev \
        libopenmpi-dev openmpi-bin \
        libproj-dev proj-bin \
        libgl1-mesa-glx \
        psfex source-extractor weightwatcher && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install uv — fast, reproducible dependency resolver and installer.
# Deps are declared in pyproject.toml; exact transitive versions are frozen
# in uv.lock. `uv sync --frozen` installs exactly what uv.lock specifies,
# so upstream changes only land when we deliberately regenerate the lockfile.
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

WORKDIR /app
COPY pyproject.toml uv.lock /app/

# Install runtime + jupyter + fitsio extras from the lockfile into /app/.venv.
# `--no-install-project` skips installing shapepipe itself (the source isn't
# copied yet); we install it `--no-deps` below once the source is available.
RUN uv sync --frozen --no-install-project --extra jupyter --extra fitsio

# Copy the source and install shapepipe into the same venv.
COPY . /app/.
RUN chown -R root:root /app && chmod -R u+rwX /app
RUN uv pip install --no-deps -e . && \
    for ext in .py .sh .bash; do \
        for script in /app/scripts/*/*$ext; do \
            link_name=$(basename $script $ext); \
            ln -s $script /usr/local/bin/$link_name; \
        done; \
    done

# Activate the uv-managed venv on container start so shapepipe_run etc
# resolve against it without explicit activation.
ENV PATH="/app/.venv/bin:${PATH}"
ENV VIRTUAL_ENV=/app/.venv
