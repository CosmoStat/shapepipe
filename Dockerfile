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
