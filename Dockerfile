# syntax=docker/dockerfile:1

# Multi-stage build for teloclip
# Stage 1: Builder - Install dependencies and build environment
FROM python:3.12-slim-bookworm AS builder

# Accept version as build argument (defaults to 0.0.0 if not provided)
ARG VERSION=0.0.0

# Set working directory
WORKDIR /build

# Install build dependencies needed for compiling Python packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    make \
    libz-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy only pyproject.toml first for better layer caching
COPY pyproject.toml .
COPY src/ ./src/
COPY README.md .
COPY LICENSE .

# Set version for setuptools-scm without requiring .git directory
ENV SETUPTOOLS_SCM_PRETEND_VERSION=${VERSION}

# Install teloclip and dependencies
# Use --no-cache-dir to reduce layer size
RUN pip install --no-cache-dir --user .

# Stage 2: Runtime - Minimal production image
FROM python:3.12-slim-bookworm AS runtime

# Add metadata labels following OCI standards
LABEL org.opencontainers.image.title="teloclip" \
    org.opencontainers.image.description="A tool for the recovery of unassembled telomeres from soft-clipped read alignments" \
    org.opencontainers.image.authors="Adam Taranto" \
    org.opencontainers.image.url="https://github.com/adamtaranto/teloclip" \
    org.opencontainers.image.source="https://github.com/adamtaranto/teloclip" \
    org.opencontainers.image.licenses="GPL-3.0-or-later"

# Copy Python packages from builder stage
COPY --from=builder /root/.local /root/.local

# Ensure scripts in .local are usable
ENV PATH=/root/.local/bin:$PATH

# Set working directory for data mounting
WORKDIR /data

# Set entrypoint to teloclip
ENTRYPOINT ["teloclip"]

# Default command shows help
CMD ["--help"]
