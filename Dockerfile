# syntax=docker/dockerfile:1
# Lumos Networks: web interface + command line, with R/Bioconductor limma.
#
# Image names are fully qualified so the same file builds with Docker and Podman.
#
#   docker build -t lumos-networks .
#   docker run --rm -p 5070:5070 -v "$PWD/lumos_jobs":/app/lumos_jobs lumos-networks     # web interface
#   docker run --rm -v "$PWD":/data -w /data lumos-networks \
#       python /app/lumos.py dataset --geo GSE117769 --outdir data                      # command line
#   docker run --rm lumos-networks python -m pytest -q /app/tests                       # tests

ARG PYTHON_IMAGE=docker.io/library/python:3.11-slim-bookworm

# ---- builder: Python dependencies, isolated from the runtime image
FROM ${PYTHON_IMAGE} AS builder
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 PIP_DISABLE_PIP_VERSION_CHECK=1
WORKDIR /app
COPY requirements.txt ./
RUN python -m pip install --no-cache-dir --upgrade pip && \
    python -m pip install --no-cache-dir --prefix=/install -r requirements.txt pytest

# ---- runtime: no build toolchain and no pip cache
FROM ${PYTHON_IMAGE}
LABEL org.opencontainers.image.title="Lumos Networks" \
      org.opencontainers.image.description="Transcriptomic analysis, functional context and prior-knowledge networks (EvoMol-Lab/BioME-UFRN)" \
      org.opencontainers.image.source="https://github.com/evomol-lab/lumos-networks" \
      org.opencontainers.image.vendor="EvoMol-Lab, BioME, UFRN" \
      org.opencontainers.image.licenses="Apache-2.0"
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONPATH=/app \
    LUMOS_HOST=0.0.0.0 \
    LUMOS_PORT=5070 \
    LUMOS_JOBS_DIR=/app/lumos_jobs
WORKDIR /app
RUN apt-get update && apt-get install -y --no-install-recommends \
        r-base-core \
        r-bioc-limma \
        r-cran-jsonlite \
        r-cran-statmod \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
COPY --from=builder /install /usr/local
COPY lumos.py lumos_core.py lumos_report.py LumosGUI.py DOCUMENTATION.md README.md LICENSE NOTICE ./
COPY r/ r/
COPY templates/ templates/
COPY static/ static/
COPY tests/ tests/
# Jobs of the web interface; mount a host directory here to keep them.
RUN useradd --create-home --uid 10001 lumos && mkdir -p lumos_jobs && chmod 777 lumos_jobs && chown -R lumos:lumos /app
USER lumos
EXPOSE 5070
STOPSIGNAL SIGINT
HEALTHCHECK --interval=30s --timeout=5s --start-period=30s --retries=3 \
  CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:5070/api/capabilities', timeout=3)" || exit 1
CMD ["python", "/app/LumosGUI.py"]
