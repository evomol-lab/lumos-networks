# syntax=docker/dockerfile:1

# ---- builder: instala as dependências Python isoladas do runtime final ----
FROM python:3.11-slim-bookworm AS builder

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /app

COPY requirements.txt ./
RUN python -m pip install --no-cache-dir --upgrade pip && \
    python -m pip install --no-cache-dir --prefix=/install -r requirements.txt

# ---- runtime: imagem final, sem toolchain de build nem cache de pip ----
FROM python:3.11-slim-bookworm

LABEL org.opencontainers.image.title="Lumos Networks" \
      org.opencontainers.image.description="Transcriptomic analysis, functional context, and prior-knowledge networks (EvoMol-Lab/BioME-UFRN)" \
      org.opencontainers.image.licenses="Apache-2.0"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    STREAMLIT_SERVER_HEADLESS=true \
    STREAMLIT_SERVER_ADDRESS=0.0.0.0 \
    STREAMLIT_SERVER_PORT=8501 \
    PYTHONPATH=/app

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
        r-base-core \
        r-bioc-limma \
        r-cran-jsonlite \
        r-cran-statmod \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

COPY --from=builder /install /usr/local

COPY . ./

RUN useradd --create-home --uid 10001 lumos && chown -R lumos:lumos /app
USER lumos

EXPOSE 8501

HEALTHCHECK --interval=30s --timeout=5s --start-period=30s --retries=3 \
  CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8501/_stcore/health', timeout=3)" || exit 1

CMD ["streamlit", "run", "Lumos_Home.py"]
