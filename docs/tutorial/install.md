# 1. Installation

## From PyPI

```bash
pip install teloclip
```

## From conda

```bash
conda install -c bioconda teloclip
```

## From source

```bash
git clone https://github.com/adamtaranto/teloclip.git
cd teloclip
pip install -e '.[dev]'
```

## With Docker

The image bundles `samtools`, because every teloclip workflow needs it.

```bash
docker pull adamtaranto/teloclip:latest

# teloclip is the entrypoint
docker run --rm -v $(pwd):/data adamtaranto/teloclip:latest --version

# override it to reach samtools
docker run --rm -v $(pwd):/data --entrypoint samtools \
  adamtaranto/teloclip:latest --version
```

See [DOCKER.md](https://github.com/adamtaranto/teloclip/blob/main/DOCKER.md) for
running whole pipelines inside one container.

## Companion tools

Teloclip does not align reads or manipulate BAM files. You will need:

- **`samtools`** — required. `extend` reads an indexed BAM and will not create
  the index itself.
- **`minimap2`** — or any long-read aligner that emits soft clips.

```bash
conda install -c bioconda samtools minimap2
```

!!! warning "Your aligner must emit soft clips, not hard clips"

    Teloclip works from the clipped sequence, and hard clips (`H` in the CIGAR)
    remove that sequence from the record entirely. Reads with hard clips are
    skipped.

    `minimap2` soft-clips primary alignments by default, which is what you want.
    Take care with `--secondary=no` plus tools that convert clips, and note that
    supplementary alignments are hard-clipped and are skipped regardless.

## Verify

```bash
teloclip --version
teloclip --help
```

You should see three sub-commands: `filter`, `extract` and `extend`.

If one is missing, teloclip will say which and why — most often a dependency
that could not be imported. `filter` and `extract` work without `pysam`;
`extend` requires it.

Next: [Mapping reads](mapping.md).
