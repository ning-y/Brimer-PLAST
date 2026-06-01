# Stage 1: Build tnBLAST (minimal build, no MPI, no NCBI toolkit)
FROM ubuntu:24.04 AS tntblast-builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    g++ \
    make \
    zlib1g-dev \
    ca-certificates \
    wget \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
RUN wget -q -O tnblast.tar.gz \
    https://github.com/jgans/thermonucleotideBLAST/archive/refs/tags/v2.77.tar.gz \
    && tar xzf tnblast.tar.gz --strip-components=1 \
    && rm tnblast.tar.gz

# Patch Makefile for minimal build (no MPI, no NCBI BLAST toolkit)
RUN sed -i Makefile \
    -e 's|^CC = mpic++|CC = g++|' \
    -e 's|-DUSE_MPI||g' \
    -e 's|-DUSE_BLAST_DB||g' \
    -e 's|^BLAST_DIR =.*|# BLAST_DIR commented out (minimal build)|'

RUN make -j$(nproc)

# Stage 2: Final image
FROM python:3.12-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    zlib1g \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=tntblast-builder /build/tntblast /usr/local/bin/tntblast

WORKDIR /app
COPY pyproject.toml /app/
COPY src/ /app/src/

RUN pip install --no-cache-dir /app

ENTRYPOINT ["brimer-plast"]
CMD ["--help"]