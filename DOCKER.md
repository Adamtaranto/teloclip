# Docker Usage Guide for Teloclip

This guide provides comprehensive instructions for using teloclip with Docker containers.

## Table of Contents

- [Quick Start](#quick-start)
- [Installation](#installation)
- [Basic Usage](#basic-usage)
- [Working with Data](#working-with-data)
- [Nextflow Integration](#nextflow-integration)
- [Building from Source](#building-from-source)
- [Troubleshooting](#troubleshooting)
- [Best Practices](#best-practices)

---

## Quick Start

Pull the latest image and run teloclip:

```bash
# Pull image
docker pull adamtaranto/teloclip:latest

# Check version
docker run --rm adamtaranto/teloclip:latest --version

# Get help
docker run --rm adamtaranto/teloclip:latest --help
```

---

## Installation

### Pull from Docker Hub

The easiest way to get teloclip is to pull the pre-built image:

```bash
docker pull adamtaranto/teloclip:latest
```

### Available Tags

- `latest` - Latest stable release
- `v0.3.2` - Specific version (recommended for reproducibility)
- `v0.3` - Latest patch in 0.3.x series
- `v0` - Latest minor in 0.x series

### Verify Installation

```bash
docker run --rm adamtaranto/teloclip:latest --version
# Output: teloclip 0.3.2
```

---

## Basic Usage

### Filter Command

Filter SAM/BAM files for telomeric overhangs:

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  adamtaranto/teloclip:latest \
  filter --ref-idx /data/ref.fa.fai /data/input.sam > output.sam
```

### Extract Command

Extract overhang sequences to FASTA files:

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  adamtaranto/teloclip:latest \
  extract --ref-idx /data/ref.fa.fai --outdir /data/overhangs /data/input.bam
```

### Extend Command

Automatically extend contigs with telomeric sequences:

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  adamtaranto/teloclip:latest \
  extend --ref /data/ref.fa --outdir /data/extended /data/overhangs.bam
```

---

## Working with Data

### Volume Mounting

Docker containers are isolated from your file system. You need to mount directories to access files:

```bash
# Mount current directory
docker run --rm -v $(pwd):/data adamtaranto/teloclip:latest filter --help

# Mount specific directory
docker run --rm -v /path/to/data:/data adamtaranto/teloclip:latest filter --ref-idx /data/ref.fa.fai /data/input.sam

# Mount multiple directories
docker run --rm \
  -v $(pwd)/input:/input:ro \
  -v $(pwd)/output:/output \
  adamtaranto/teloclip:latest \
  filter --ref-idx /input/ref.fa.fai /input/data.sam > /output/filtered.sam
```

### Using `:ro` for Read-Only Mounts

Protect input data from accidental modification:

```bash
docker run --rm \
  -v $(pwd)/input:/input:ro \
  -v $(pwd)/output:/output \
  adamtaranto/teloclip:latest \
  filter --ref-idx /input/ref.fa.fai /input/data.sam > /output/filtered.sam
```

### Working Directory

The container's working directory is `/data`. Files without absolute paths are relative to this:

```bash
# Mount data directory
docker run --rm -v $(pwd)/data:/data adamtaranto/teloclip:latest \
  filter --ref-idx ref.fa.fai input.sam > output.sam
```

---

## Pipeline Integration

### With Piping and Other Containers

Combine teloclip with other bioinformatics containers:

```bash
# Align with minimap2, filter with teloclip, sort with samtools
docker run --rm -v $(pwd):/data minimap2/minimap2 \
  -ax map-pb /data/ref.fa /data/reads.fq.gz | \
docker run --rm -i -v $(pwd):/data adamtaranto/teloclip:latest \
  filter --ref-idx /data/ref.fa.fai /dev/stdin | \
docker run --rm -i -v $(pwd):/data biocontainers/samtools:v1.18 \
  sort > /data/overhangs.bam
```

### Shell Script Example

Create a complete workflow script:

```bash
#!/usr/bin/env bash
# teloclip-pipeline.sh

set -euo pipefail

REF="ref.fa"
READS="reads.fq.gz"
OUTDIR="results"
MOTIF="TTAGGG"

mkdir -p ${OUTDIR}

# Index reference
docker run --rm -v $(pwd):/data \
  biocontainers/samtools:v1.18 \
  faidx /data/${REF}

# Align reads
docker run --rm -v $(pwd):/data \
  minimap2/minimap2 -ax map-pb /data/${REF} /data/${READS} \
  > ${OUTDIR}/aligned.sam

# Filter for telomeric overhangs
docker run --rm -v $(pwd):/data \
  adamtaranto/teloclip:latest \
  filter --ref-idx /data/${REF}.fai \
  --motifs ${MOTIF} \
  /data/${OUTDIR}/aligned.sam \
  > ${OUTDIR}/overhangs.sam

# Convert to BAM
docker run --rm -v $(pwd):/data \
  biocontainers/samtools:v1.18 \
  sort /data/${OUTDIR}/overhangs.sam \
  > ${OUTDIR}/overhangs.bam

# Extract sequences
docker run --rm -v $(pwd):/data \
  adamtaranto/teloclip:latest \
  extract --ref-idx /data/${REF}.fai \
  --outdir /data/${OUTDIR}/extracted \
  /data/${OUTDIR}/overhangs.bam

echo "Pipeline complete! Results in ${OUTDIR}/"
```

---

## Nextflow Integration

See [examples/nextflow/](../examples/nextflow/) for complete Nextflow workflows and modules.

### Quick Example

```groovy
process TELOCLIP_FILTER {
    container 'adamtaranto/teloclip:latest'

    input:
    path sam
    path fai

    output:
    path "overhangs.sam"

    script:
    """
    teloclip filter --ref-idx ${fai} ${sam} > overhangs.sam
    """
}
```

### Multi-Container Best Practice

```groovy
// Each tool in separate container
process BAM_TO_SAM {
    container 'biocontainers/samtools:v1.18'
    input: path bam
    output: path "input.sam"
    script: "samtools view -h ${bam} > input.sam"
}

process TELOCLIP_FILTER {
    container 'adamtaranto/teloclip:latest'
    input: path sam, path fai
    output: path "filtered.sam"
    script: "teloclip filter --ref-idx ${fai} ${sam} > filtered.sam"
}
```

---

## Building from Source

### Local Build

```bash
# Clone repository
git clone https://github.com/adamtaranto/teloclip.git
cd teloclip

# Build image
docker build -t teloclip:local .

# Test build
docker run --rm teloclip:local --version
```

### Multi-Architecture Build

```bash
# Use provided script
./scripts/build-docker.sh

# Or manually with buildx
docker buildx build \
  --platform linux/amd64,linux/arm64 \
  --tag teloclip:multi \
  --load \
  .
```

### Build Options

```bash
# Build specific version
docker build -t teloclip:0.3.2 --build-arg VERSION=0.3.2 .

# Build with cache
docker build --cache-from teloclip:latest -t teloclip:dev .
```

---

## Troubleshooting

### Common Issues

#### 1. Permission Denied

**Problem**: Cannot write to output directory

```text
docker: Error response from daemon: error while creating mount source path '/data/output': permission denied
```

**Solution**: Check directory permissions or use absolute paths

```bash
# Create directory first
mkdir -p output
chmod 755 output

# Or use absolute path
docker run --rm -v $(pwd)/output:/output adamtaranto/teloclip:latest ...
```

#### 2. File Not Found

**Problem**: Cannot find input files

```text
FileNotFoundError: [Errno 2] No such file or directory: '/data/input.sam'
```

**Solution**: Ensure files are in mounted directory

```bash
# Check file exists
ls -l data/input.sam

# Mount parent directory
docker run --rm -v $(pwd)/data:/data adamtaranto/teloclip:latest \
  filter --ref-idx /data/ref.fa.fai /data/input.sam
```

#### 3. Image Not Found

**Problem**: Cannot pull or run image

```text
Unable to find image 'adamtaranto/teloclip:latest' locally
docker: Error response from daemon: pull access denied
```

**Solution**: Check image name and Docker Hub access

```bash
# Verify image name
docker search adamtaranto/teloclip

# Pull explicitly
docker pull adamtaranto/teloclip:latest

# Check Docker daemon
docker version
```

#### 4. Out of Memory

**Problem**: Container killed due to memory

```text
Killed
```

**Solution**: Increase Docker memory limit

```bash
# Docker Desktop: Preferences → Resources → Memory
# Or run with memory limit
docker run --rm --memory="8g" -v $(pwd):/data \
  adamtaranto/teloclip:latest filter ...
```

#### 5. Platform Mismatch

**Problem**: Running AMD64 image on ARM64 (or vice versa)

```text
WARNING: The requested image's platform (linux/amd64) does not match the detected host platform
```

**Solution**: Multi-arch images should work via emulation, but specify if needed

```bash
# Force specific platform
docker run --rm --platform linux/amd64 adamtaranto/teloclip:latest --version
```

---

## Best Practices

### 1. Pin Versions for Reproducibility

```bash
# Good - specific version
docker pull adamtaranto/teloclip:v0.3.2

# Avoid - may change
docker pull adamtaranto/teloclip:latest
```

### 2. Use Read-Only Mounts for Input

```bash
docker run --rm \
  -v $(pwd)/input:/input:ro \
  -v $(pwd)/output:/output \
  adamtaranto/teloclip:latest filter --ref-idx /input/ref.fa.fai /input/data.sam
```

### 3. Clean Up Containers

```bash
# Always use --rm for automatic cleanup
docker run --rm adamtaranto/teloclip:latest --help

# Manually clean if needed
docker container prune
```

### 4. Use Absolute Paths

```bash
# Good
docker run --rm -v /full/path/to/data:/data adamtaranto/teloclip:latest ...

# Also good
docker run --rm -v $(pwd)/data:/data adamtaranto/teloclip:latest ...
```

### 5. Check Exit Codes

```bash
#!/bin/bash
set -e  # Exit on error

docker run --rm -v $(pwd):/data adamtaranto/teloclip:latest \
  filter --ref-idx /data/ref.fa.fai /data/input.sam > output.sam

if [ $? -eq 0 ]; then
  echo "Success!"
else
  echo "Failed!" >&2
  exit 1
fi
```

### 6. Log Output

```bash
# Redirect stdout and stderr
docker run --rm -v $(pwd):/data adamtaranto/teloclip:latest \
  filter --ref-idx /data/ref.fa.fai /data/input.sam \
  > output.sam 2> teloclip.log
```

### 7. Use Environment Variables

```bash
# Set common variables
export DATA_DIR=/path/to/data
export TELOCLIP_IMAGE=adamtaranto/teloclip:v0.3.2

docker run --rm -v ${DATA_DIR}:/data ${TELOCLIP_IMAGE} \
  filter --ref-idx /data/ref.fa.fai /data/input.sam
```

---

## Performance Considerations

### Docker vs Native

- **Docker**: ~2-5% overhead for I/O
- **Native**: Fastest, but requires Python environment
- **Recommendation**: Docker for reproducibility, native for HPC

### Volume Performance

```bash
# Fastest: Native Docker volumes
docker volume create teloclip-data
docker run --rm -v teloclip-data:/data adamtaranto/teloclip:latest ...

# Good: Bind mounts with delegated mode (macOS)
docker run --rm -v $(pwd):/data:delegated adamtaranto/teloclip:latest ...

# Slower: Default bind mounts
docker run --rm -v $(pwd):/data adamtaranto/teloclip:latest ...
```

### Resource Limits

```bash
# Set memory limit
docker run --rm --memory="4g" ...

# Set CPU limit
docker run --rm --cpus="2" ...

# Both
docker run --rm --memory="4g" --cpus="2" ...
```

---

## Security

### Run as Non-Root (Advanced)

```dockerfile
# Custom Dockerfile
FROM adamtaranto/teloclip:latest
USER 1000:1000
```

### Scan for Vulnerabilities

```bash
# Using Docker scan
docker scan adamtaranto/teloclip:latest

# Using Trivy
trivy image adamtaranto/teloclip:latest
```

---

## Additional Resources

- **GitHub Repository**: <https://github.com/adamtaranto/teloclip>
- **Docker Hub**: <https://hub.docker.com/r/adamtaranto/teloclip>
- **Issues**: <https://github.com/adamtaranto/teloclip/issues>
- **Nextflow Examples**: [examples/nextflow/](../examples/nextflow/)

---

## Getting Help

If you encounter issues:

1. Check this troubleshooting guide
2. Search existing [GitHub issues](https://github.com/adamtaranto/teloclip/issues)
3. Create a new issue with:
   - Docker version (`docker --version`)
   - Image tag used
   - Complete command run
   - Error message
   - Operating system

---

**Last Updated**: November 6, 2025
**Docker Image Version**: 0.3.2
