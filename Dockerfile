# syntax=docker/dockerfile:1
# Lumos Networks: web interface + command line, with R/Bioconductor limma.
#
# Image names are fully qualified so the same file builds with Docker and Podman.
# Podman needs --format docker to keep the HEALTHCHECK (its default OCI format drops it).
#
#   podman build --format docker -t docker.io/evomol/lumos-networks:2.0.0 .
#   podman run --rm -p 5070:5070 docker.io/evomol/lumos-networks:2.0.0                  # web interface
#   podman run --rm -v "$PWD":/data -w /data docker.io/evomol/lumos-networks:2.0.0 \
#       python /opt/lumos/lumos.py dataset --geo GSE117769 --outdir data            # command line
#   podman run --rm docker.io/evomol/lumos-networks:2.0.0 \
#       python -m pytest -q /opt/lumos/tests                                        # tests

ARG PYTHON_IMAGE=docker.io/library/python:3.11-slim-bookworm

# ---- build: Python dependencies in a virtual environment, without pip caches
FROM ${PYTHON_IMAGE} AS build
ENV PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1
RUN python -m venv /opt/venv
COPY requirements.txt /tmp/requirements.txt
RUN /opt/venv/bin/pip install --upgrade pip && \
    /opt/venv/bin/pip install -r /tmp/requirements.txt pytest

# ---- runtime
FROM ${PYTHON_IMAGE}
ARG VERSION=2.0.0
LABEL org.opencontainers.image.title="Lumos Networks" \
      org.opencontainers.image.version="${VERSION}" \
      org.opencontainers.image.description="Transcriptomic analysis, functional context and prior-knowledge networks: web interface and command line, with limma (R) and PyDESeq2" \
      org.opencontainers.image.url="https://hub.docker.com/r/evomol/lumos-networks" \
      org.opencontainers.image.source="https://github.com/evomol-lab/lumos-networks" \
      org.opencontainers.image.documentation="https://github.com/evomol-lab/lumos-networks/blob/main/DOCUMENTATION.md" \
      org.opencontainers.image.vendor="EvoMol-Lab, BioME, UFRN" \
      org.opencontainers.image.licenses="Apache-2.0"
ENV PATH=/opt/venv/bin:$PATH \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPATH=/opt/lumos \
    MPLBACKEND=Agg \
    LUMOS_HOST=0.0.0.0 \
    LUMOS_PORT=5070 \
    LUMOS_JOBS_DIR=/opt/lumos/lumos_jobs
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        r-base-core \
        r-bioc-limma \
        r-cran-jsonlite \
        r-cran-statmod \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
COPY --from=build /opt/venv /opt/venv
WORKDIR /opt/lumos
COPY lumos.py lumos_core.py lumos_report.py LumosGUI.py DOCUMENTATION.md README.md LICENSE NOTICE ./
COPY r/ r/
COPY templates/ templates/
COPY static/ static/
COPY tests/ tests/
# Jobs of the web interface; mount a host directory here to keep them. The image runs as root,
# like the other EvoMol images, so a mounted folder is writable with rootless Podman; with rootful
# Docker, add --user "$(id -u):$(id -g)" so the files belong to you.
RUN mkdir -p lumos_jobs && chmod 777 lumos_jobs
EXPOSE 5070
# Python running as PID 1 ignores SIGTERM; SIGINT stops the Flask server cleanly.
STOPSIGNAL SIGINT
HEALTHCHECK --interval=30s --timeout=5s --start-period=30s --retries=3 \
  CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:5070/api/capabilities', timeout=3)" || exit 1
# Absolute path: Apptainer images built from this one (docker://...) ignore WORKDIR.
CMD ["python", "/opt/lumos/LumosGUI.py"]
