# syntax=docker/dockerfile:1.7
#
# Two-target image:
#   --target runtime  →  minimal, for canfar batch jobs and downstream stacks
#   --target dev      →  runtime + everyday CLI tools + all extras (test,
#                        lint, doc, …); default if --target is omitted
#
# Both share the `base` stage (system deps + uv + lockfile copy), so the
# heavy apt + wheel-resolution work is cached once.

# ----------------------------------------------------------------------
# base — system deps shared by every target
# ----------------------------------------------------------------------
FROM python:3.12-slim-bookworm AS base

LABEL maintainer="martin.kilbinger@cea.fr"

ENV SHELL=/bin/bash \
    QT_QPA_PLATFORM=offscreen \
    PIP_NO_CACHE_DIR=1 \
    DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    # Make the image well-behaved on read-only filesystems (apptainer/SIF):
    #   UV_NO_SYNC      — `uv run` skips the auto-sync that would otherwise
    #                     mutate /app/.venv. Run with what's baked in.
    #   UV_CACHE_DIR    — uv's package cache. /tmp is tmpfs, writable in
    #                     read-only SIFs and writable sandboxes alike.
    #   COVERAGE_FILE   — pytest-cov writes its data file here instead of
    #                     /app/.coverage; needed because pyproject sets
    #                     `addopts = --cov=shapepipe` and pytest-cov erases
    #                     the file at startup.
    UV_NO_SYNC=1 \
    UV_CACHE_DIR=/tmp/uv-cache \
    COVERAGE_FILE=/tmp/.coverage

# System dependencies — three categories:
#  - astromatic binaries (psfex, source-extractor, weightwatcher) ship as
#    Debian packages on bookworm; preferred over building from source.
#  - compilers and dev libs needed to build the heavier wheels (galsim,
#    mpi4py, python-pysap, fitsio).
#  - libgl1, proj, fftw at runtime for skyproj/PyQt5/galsim.
# OpenMPI is deliberately NOT installed from Debian here — bookworm ships
# OpenMPI 4.1.4 / PMIx 2.x, which breaks hybrid MPI on modern clusters. It is
# built from source in the next stanza; see there for the full reasoning.
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
        libproj-dev proj-bin \
        libgl1-mesa-glx \
        psfex source-extractor weightwatcher && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# OpenMPI from source — required for hybrid Apptainer MPI on HPC clusters.
#
# On a cluster ShapePipe runs as a standard Apptainer "hybrid" MPI job: the
# host's `mpirun` launches one container rank per slot, and the OpenMPI inside
# the image wires the ranks together through PMIx. That handshake requires the
# container's PMIx to be compatible with the host launcher's. Debian bookworm's
# package is OpenMPI 4.1.4 with PMIx 2.x; modern clusters (e.g. candide) now run
# OpenMPI 5.0.x with PMIx 5.x, and a PMIx 2 client cannot talk to a PMIx 5
# server — so every rank silently degrades to a standalone "rank 0 of 1" and the
# job runs N independent copies instead of one N-rank job. Building OpenMPI
# 5.0.x here (with its bundled PMIx 5 / PRRTE) matches those hosts; the 5.0.x
# series is mutually PMIx-compatible, so this image works against any host
# openmpi/5.0.x module. The stock mpi4py wheel (from uv.lock) dlopens
# libmpi.so.40, the soname this build provides, so it needs no rebuild.
#
# --disable-dlopen links every MCA component statically into libmpi / libpmix:
# it sidesteps an internal-openpmix configure failure (the `pdl` component wants
# libltdl headers otherwise) and is the right posture for a container anyway —
# no fragile runtime dlopen of plugin .so files across the SIF / bind boundary.
ARG OMPI_VERSION=5.0.8
ARG OMPI_SERIES=v5.0
RUN cd /tmp && \
    wget -q "https://download.open-mpi.org/release/open-mpi/${OMPI_SERIES}/openmpi-${OMPI_VERSION}.tar.bz2" && \
    tar xjf "openmpi-${OMPI_VERSION}.tar.bz2" && \
    cd "openmpi-${OMPI_VERSION}" && \
    ./configure --prefix=/opt/ompi \
        --with-pmix=internal --with-prrte=internal \
        --with-hwloc=internal --with-libevent=internal \
        --disable-dlopen --disable-sphinx && \
    make -j"$(nproc)" && make install && \
    cd / && rm -rf /tmp/openmpi-*
ENV PATH="/opt/ompi/bin:${PATH}" \
    LD_LIBRARY_PATH="/opt/ompi/lib:${LD_LIBRARY_PATH}"

# uv — fast reproducible Python deps installer. pyproject.toml + uv.lock
# are the SSOT; `uv sync --frozen` installs exactly what uv.lock specifies,
# so upstream changes only land when we deliberately regenerate the lockfile.
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

WORKDIR /app
COPY pyproject.toml uv.lock /app/

# ----------------------------------------------------------------------
# runtime — minimal target for batch jobs and downstream FROM clauses
# ----------------------------------------------------------------------
FROM base AS runtime
LABEL description="ShapePipe runtime — slim Python + uv-frozen deps"

# Lockfile-frozen Python deps + jupyter + fitsio. Test/lint/doc extras
# are intentionally left out here; they live in the dev target.
RUN uv sync --frozen --no-install-project --extra jupyter --extra fitsio

# Copy the source and install shapepipe into the same venv.
COPY . /app/.
# go+rwX so non-root users on canfar/skaha can read/traverse /app and
# write into the venv when they need to (e.g. uv add for ad-hoc deps).
RUN chmod -R go+rwX /app && \
    uv pip install --no-deps -e . && \
    for ext in .py .sh .bash; do \
        for script in /app/scripts/*/*$ext; do \
            link_name=$(basename $script $ext); \
            ln -s $script /usr/local/bin/$link_name; \
        done; \
    done

# Activate the uv-managed venv on container start so shapepipe_run etc
# resolve against it without explicit activation.
ENV PATH="/app/.venv/bin:${PATH}" \
    VIRTUAL_ENV=/app/.venv

# ----------------------------------------------------------------------
# dev — everyday working environment (default target)
# ----------------------------------------------------------------------
FROM base AS dev
LABEL description="ShapePipe dev — runtime + interactive CLI tools + all extras"

# Interactive tools for working inside the container. Curated subset of
# cailmdaley/containers focused on the search/edit/process loop; heavier
# tooling (neovim, polspice, quarto, zellij) is intentionally not here.
RUN apt-get update -y --quiet && \
    apt-get install -y --no-install-recommends \
        vim \
        less \
        tmux \
        htop \
        procps \
        ripgrep \
        fd-find \
        jq \
        bat \
        curl \
        ca-certificates \
        git-lfs \
        rsync \
        unzip zip \
        openssh-client \
        locales && \
    if command -v batcat >/dev/null; then ln -sf "$(command -v batcat)" /usr/local/bin/bat; fi && \
    if command -v fdfind >/dev/null; then ln -sf "$(command -v fdfind)" /usr/local/bin/fd; fi && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# All extras pre-installed (dev = doc + jupyter + lint + release + test +
# fitsio). Pre-installing avoids the read-only-fs failure Martin hit when
# trying to live `uv sync --extra test` inside the runtime image on canfar.
RUN uv sync --frozen --no-install-project --extra dev

COPY . /app/.
RUN chmod -R go+rwX /app && \
    uv pip install --no-deps -e . && \
    for ext in .py .sh .bash; do \
        for script in /app/scripts/*/*$ext; do \
            link_name=$(basename $script $ext); \
            ln -s $script /usr/local/bin/$link_name; \
        done; \
    done

ENV PATH="/app/.venv/bin:${PATH}" \
    VIRTUAL_ENV=/app/.venv
