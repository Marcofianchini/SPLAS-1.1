FROM python:3.12-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends gfortran make libblas-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/splas
COPY . .
RUN python -m pip install --no-cache-dir numpy scipy matplotlib meson ninja setuptools wheel \
    && make build \
    && python -m pip install --no-cache-dir -e .

ENTRYPOINT ["splas"]
CMD ["--help"]
