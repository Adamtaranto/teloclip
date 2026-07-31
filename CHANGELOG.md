# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

- **Contig shortening**: `teloclip extend` could silently make a contig *shorter*. Extending an
  end trims the bases the supporting read did not cover before grafting on its soft clip, but
  nothing checked that the clip reached past the contig terminus. A read aligning from position
  10 with a 3 bp clip trimmed 9 bases and added 3. `validate_extension` did not catch it because
  it was handed the already-trimmed sequence as the "original", making its length check trivially
  true. Extensions are now rejected unless they add at least one net base.
- **`teloclip extract` crash**: `EnhancedStreamingSamFilter` read `self.exclude_secondary`
  without ever assigning it, raising `AttributeError` on the first secondary alignment in the
  input. Only inputs already stripped of secondary alignments ran to completion.
- **`teloclip extract` motif filtering**: a read whose left clip contained no motif skipped the
  whole record, so its right clip was never examined. Each end is now judged independently.

### Changed

Each sub-command carried its own copy of the "is this soft clip a valid terminal overhang?"
test, and the three copies disagreed. They are now a single shared implementation
(`teloclip.overhang`) using 1-based inclusive coordinates. This changes which reads are
accepted:

- **`alignment_end` is now inclusive everywhere.** `filter` and `extract` computed it as
  `POS + reference_span`, one base past the last aligned base. Their right-end tolerance was
  therefore one base wider than intended, and out of step with their own left end.
- **`max_break` is now an inclusive gap, applied identically at both ends.** `gap <= max_break`,
  where the gap is the count of terminal contig bases the alignment does not cover.
  `--max-break 0` now means the alignment must reach the terminal base exactly.
- **`min_clip` now measures bases past the terminus, not raw clip length**, matching what its
  help text has always said. It was previously applied to the raw clip on the left end and
  **not enforced at all on the right end**, where the test reduced to "overhang >= 0" and so
  accepted clips that reached the terminus without passing it.
- **Anchor length in `extend` now counts `M`/`=`/`X` only.** It was computed as read length
  minus clips, which counts insertions as anchoring evidence. `filter` and `extract` were
  already correct. Expect a small tightening on insertion-rich long-read alignments.
- **`teloclip extend --max-break` default moves from 10 to 50**, matching `filter` and
  `extract`. Users previously filtered at 50 and then silently lost candidates at 10.
- **Overhang selection now ranks on net gain, preferring less trimming.** A long clip anchored
  well inside a contig can contribute less novel sequence than a shorter clip flush with the
  terminus. Where two candidates offer comparable gain, the one that discards less assembly
  wins, so a single extra base cannot justify replacing polished consensus with one raw read.
  `--min-extension` likewise now gates on novel bases contributed.
- **Dry runs now account for trimming.** They hard-coded the trim to zero, so a dry run
  predicted a longer contig than the real run produced wherever trimming occurred.

- **Output contig order now matches input order.** `extend` wrote extended contigs first and
  appended everything else, so excluded, unsupported and failed contigs were all shunted to the
  end of the FASTA.
- **Extension amounts are now reported net of trimming.** `Original bp + Total +bp = Final bp`
  is now an identity you can check; previously the report showed raw clip lengths and did not
  reconcile with the sequences it described.
- **Contigs with anomalous overhang coverage are reported, not excluded.** `--exclude-outliers`
  silently dropped contigs and the report was passed empty outlier lists regardless, so the
  exclusion was invisible. Detection also used mean and standard deviation (both inflated by
  the outliers being sought), flagged both tails, and could not fire at all on small
  assemblies. It now uses median/MAD on the high tail with a minimum-contig floor, and only
  advises. `--exclude-outliers` is accepted but ignored, with a warning; `--outlier-threshold`
  now sets flagging sensitivity and defaults to 3.5.
- **Unknown `--exclude-contigs` names are now an error.** A misspelled name previously matched
  nothing and the run continued, extending a contig the user believed was held back.
- **`filter` exclusion counts now reconcile.** Reads with no usable soft clip fell through every
  bucket while still counting toward the discard total. The two near-identical summary blocks
  are merged into one.
- **`--match-anywhere` no longer applies `--min-repeats` twice**, which demanded the repeat
  count squared.

### Added

- **`teloclip extend --min-clip`** (default 1, clamped to at least 1), so the
  `(max_break, min_clip, min_anchor)` triple is uniform across all three sub-commands.
- **`--logfile`** on all three sub-commands. Writing logs to a file was fully implemented but
  no CLI exposed it.
- **`teloclip extend --overhang-log`**, a TSV of every accepted overhang: contig, end,
  alignment span, gap from the terminus, clip length, overhang length and anchor length.
- **Per-contig overhang counts and a distribution histogram** logged by `filter` and `extend`.
- **`teloclip extend --html-report`**, a self-contained HTML report showing the evidence
  behind each extension: every supporting read laid out against the contig terminus it
  supports, with anchor and soft clip distinguished, motif matches highlighted, the selected
  read marked, and read details on hover. Plus a strip plot of overhang depth across contigs
  with a median reference line and flagged ends ringed and labelled, backed by a table view.
  `--html-max-reads` (default 25) caps the rows rendered per contig end.

  Reads are placed by walking their CIGAR, so indels shift them against the
  reference as the aligner recorded: a deletion gaps the read where the contig
  has bases, and an insertion in one read opens a column in every other row
  including the contig, which is what keeps the rows aligned. The anchored
  region shown is at least 500 bases, or the longest anchor at that end if
  longer. A position ruler runs in contig bases with 0 at the terminus, reads
  overhanging both ends of a contig are marked, hovering a read gives a table of
  its SAM record, and chart points show their contig on hover or click.

### Removed

- **`analysis.collect_overhang_stats`**, an unused fourth copy of the overhang test carrying
  its own hard-coded thresholds.
- **`streaming_io.copy_unmodified_contigs`**, replaced by a single ordered write pass.

---

## [0.3.4] - 2025-11-07

### Added

- **Docker containerization**: Production-ready multi-stage Dockerfile with multi-architecture support (amd64/arm64)
- **Docker Hub integration**: Automated CI/CD pipeline via GitHub Actions for building and publishing images
- **Nextflow workflow**: Complete example pipeline (`examples/nextflow/teloclip.nf`) with nf-core style modules for filter, extract, and extend
- **Docker documentation**: Comprehensive `DOCKER.md` guide with usage examples and best practices
- **Build scripts**: `scripts/build-docker.sh` and `scripts/test-docker.sh` for local development and testing

### Changed

- **License**: Migrated from MIT to GPL-3.0-or-later for better alignment with open science principles
- **Docker workflow**: GitHub Actions now triggers only on version tags and pull requests (not on regular commits to main)
- **Version management**: Implemented PEP 440 compliant version conversion in build scripts

### Fixed

- **Nextflow modules**: Corrected `teloclip extract` to accept SAM input (not BAM)
- **Nextflow modules**: Fixed version parsing in all modules to match actual CLI output format
- **Nextflow workflow**: Fixed SAM_TO_BAM process to properly output separate SAM file for extract step

---

## [0.3.3]

### 🐳 Docker & Container Infrastructure

#### Docker Support

- **Production-ready Dockerfile**: Multi-stage build using `python:3.12-slim-bookworm` base image
  - Optimized build process with separate builder and runtime stages
  - Minimal final image size (~80-100MB) with only production dependencies
  - Smart layer caching for faster rebuilds
- **Multi-architecture support**: Images built for both `linux/amd64` and `linux/arm64` platforms
- **Automated version injection**: PEP 440 compliant versions extracted from git tags and passed as build arguments
- **Build optimization**: Comprehensive `.dockerignore` excluding tests, docs, and development files

#### Docker Hub Integration

- **Automated publishing**: GitHub Actions workflow for building and publishing Docker images
- **Smart tagging strategy**:
  - `latest` - Most recent release from main branch
  - `v{major}.{minor}.{patch}` - Full semantic version tags
  - `v{major}.{minor}` - Minor version tags
  - `v{major}` - Major version tags
  - `pr-{number}` - Pull request preview builds
- **CI/CD pipeline**: Automated builds on push to main, tags, and pull requests
- **PR testing**: Automatic image validation for pull requests with `--version` and `--help` tests
- **Build caching**: GitHub Actions cache for faster subsequent builds

#### Development Tools

- **Local build script** (`scripts/build-docker.sh`): Helper script for local multi-architecture builds
  - Automatic version detection from git tags
  - Configurable platforms and image tags
  - Support for push to registry or local load
- **Test automation** (`scripts/test-docker.sh`): Comprehensive Docker image validation
  - Version and help command tests
  - Subcommand functionality validation (filter/extract/extend)
  - Image size verification
  - Test data processing validation

### 🔄 Nextflow Pipeline Integration

#### Complete Workflow Examples

- **Example pipeline** (`examples/nextflow/teloclip.nf`): End-to-end telomere analysis workflow
  - Reference indexing with samtools
  - BAM/SAM conversion and filtering
  - Teloclip filtering, extraction, and extension steps
  - Proper container isolation (each tool in its own container)
  - Comprehensive parameter handling

#### nf-core Style Modules

- **TELOCLIP_FILTER module**: Filter SAM/BAM files for terminal soft-clipped alignments
  - Configurable motif matching with `ext.args`
  - Meta map support for sample tracking
  - Version tracking output
- **TELOCLIP_EXTRACT module**: Extract overhang sequences to FASTA files
  - Per-contig organization with prefix support
  - Optional statistics and motif counting
  - Buffered output for large datasets
- **TELOCLIP_EXTEND module**: Automatically extend contigs with telomeric sequences
  - Takes filtered BAM + reference genome as input
  - Outputs extended genome FASTA and statistics report
  - Supports dry-run mode and outlier detection

#### Pipeline Documentation

- **Complete integration guide** (`examples/nextflow/README.md`):
  - Quick start examples
  - Module usage patterns
  - Container strategy explanation
  - Parameter reference
  - Output structure documentation
  - Troubleshooting tips
  - nf-core integration guidelines

### 📚 Documentation Enhancements

#### Docker Documentation

- **Comprehensive Docker guide** (`DOCKER.md`): 557-line complete reference
  - Quick start instructions
  - Installation methods for all tag types
  - Usage examples for all subcommands (filter/extract/extend)
  - Volume mounting patterns
  - Pipeline integration with shell scripts
  - Nextflow workflow examples
  - Building from source instructions
  - Common troubleshooting scenarios
  - Best practices and security considerations
  - Performance optimization tips

#### README Updates

- **Docker badges**: Added Docker Hub version and pull count badges
- **Docker installation section**: New installation method #5 with quick start commands
- **Cross-references**: Links to DOCKER.md and Nextflow examples
- **Container-first approach**: Documentation now emphasizes Docker as a primary installation method

### ⚖️ License Change

#### Migration to GPL-3.0-or-later

- **License update**: Changed from MIT to GPL-3.0-or-later
- **Rationale**: Better alignment with open science principles and academic software sharing
- **Compatibility**: Maintains compatibility with derivative works and academic use
- **Documentation**: Updated all references in:
  - `LICENSE` file
  - `README.md` badges
  - `pyproject.toml` metadata
  - `Dockerfile` OCI labels
  - Nextflow module metadata

### 🔧 Infrastructure Improvements

#### Build System Enhancements

- **Version management**: PEP 440 compliant version conversion in build scripts
  - Converts git describe output (`v0.3.2-6-gd5ce6fc-dirty`) to valid Python versions (`0.3.2.post6+gd5ce6fc.dirty`)
  - Supports release tags, post-releases, and local versions
- **Environment variables**: `SETUPTOOLS_SCM_PRETEND_VERSION` for builds without .git directory
- **Build reproducibility**: Version baked into Docker images at build time

#### GitHub Actions

- **Docker build workflow** (`.github/workflows/docker-build.yml`):
  - QEMU setup for multi-architecture builds
  - Docker Buildx configuration
  - Docker Hub authentication with secrets
  - Metadata extraction with semantic versioning
  - Conditional push/load based on event type
  - Job summaries with build information
  - PR-specific single-platform builds for testing

### 🎯 Container Best Practices

#### Security & Optimization

- **Minimal attack surface**: Runtime image contains only essential Python packages
- **No development tools**: Build tools (gcc, g++, make) excluded from final image
- **Reproducible builds**: Locked Python version and base image
- **OCI labels**: Complete metadata following Open Container Initiative standards
- **Non-root execution**: Container runs with appropriate user permissions

#### Workflow Integration

- **Single-tool containers**: Each tool (teloclip, samtools) in separate containers
- **Follows nf-core best practices**: Module structure compatible with nf-core pipelines
- **Version pinning**: Explicit container versions for reproducibility
- **Data mounting**: Clear patterns for volume mounting in examples

### 📦 Repository Structure

#### New Files

```text
.dockerignore                              # Docker build optimization
Dockerfile                                 # Multi-stage production build
DOCKER.md                                  # Complete Docker documentation
scripts/
  ├── build-docker.sh                     # Multi-arch build helper
  └── test-docker.sh                      # Container validation tests
.github/workflows/
  └── docker-build.yml                    # Automated Docker CI/CD
examples/nextflow/
  ├── teloclip.nf                         # Complete example workflow
  ├── README.md                           # Integration documentation
  └── modules/teloclip/
      ├── filter/
      │   ├── main.nf                     # Filter module
      │   └── meta.yml                    # Module metadata
      ├── extract/
      │   ├── main.nf                     # Extract module
      │   └── meta.yml                    # Module metadata
      └── extend/
          ├── main.nf                     # Extend module
          └── meta.yml                    # Module metadata
```

### 🚀 Usage Examples

#### Docker Quick Start

```bash
# Pull latest image
docker pull adamtaranto/teloclip:latest

# Check version
docker run --rm adamtaranto/teloclip:latest --version

# Filter alignments with volume mounting
docker run --rm -v $(pwd):/data adamtaranto/teloclip:latest filter \
  --ref-idx /data/ref.fa.fai /data/input.sam > filtered.sam
```

#### Nextflow Integration

```bash
# Run complete workflow
nextflow run examples/nextflow/teloclip.nf \
  --bam input.bam \
  --ref reference.fa \
  --motifs TTAGGG \
  --outdir results
```

### 🔄 Breaking Changes

None - This release is fully backward compatible with existing CLI usage.

### 📝 Migration Notes

#### For Docker Users

- Docker images are now the recommended installation method for production workflows
- No changes required for existing users - pip and conda installation still fully supported
- Docker provides better reproducibility and easier dependency management

#### For Pipeline Developers

- Nextflow modules follow nf-core conventions and can be integrated into existing pipelines
- Container strategy uses separate containers per tool (best practice)
- All modules include comprehensive metadata for pipeline integration

---

## [0.3.2] - 2025-10-26

### 🚀 Major Features Added

#### New Automatic Contig Extension (`teloclip extend`)

- **Intelligent genome completion**: Automatically extend draft genome contigs using telomeric overhang sequences from soft-clipped alignments
- **Quality control system**: Built-in outlier detection to prevent extension of circular genomes (mitochondria, chloroplasts)
- **Dry-run capability**: Preview proposed extensions before applying changes with `--dry-run` flag
- **Comprehensive validation**: Multi-level validation including homopolymer run detection and anchor quality assessment
- **Flexible exclusion**: Skip specific contigs by name (`--exclude-contig`) or from file lists (`--exclude-contig-file`)

#### Complete CLI Architecture Redesign

- **Modular sub-commands**: Replaced single-purpose CLI with three specialized commands:
  - `teloclip filter` - Enhanced alignment filtering with rich logging
  - `teloclip extract` - Advanced sequence extraction with multiple output formats
  - `teloclip extend` - **NEW** automatic contig extension functionality
- **Unified interface**: Consistent argument patterns and help documentation across all commands
- **Backward compatibility**: Maintains core functionality while modernizing interface

#### Memory-Efficient Streaming Architecture

- **Large genome support**: Process multi-gigabase genomes without loading entire files into memory
- **Indexed file access**: Efficient BAM/FASTA processing using pysam for contig-by-contig operations
- **Buffered I/O operations**: Optimized file writing with configurable buffer sizes

### ✨ Enhanced Features

#### Advanced Motif Analysis

- **Fuzzy matching support**: Configurable tolerance for motif variations (±1 character in homopolymer runs)
- **Comprehensive pattern counting**: Detailed motif occurrence statistics in overhang regions
- **Optimized regex processing**: Pre-compiled patterns for improved performance
- **Automatic reverse complement**: Seamless handling of both forward and reverse motif orientations

#### Rich Statistical Analysis and Reporting

- **Overhang statistics**: Detailed analysis with Z-score based outlier detection
- **Extension reports**: Human-readable summaries of proposed and applied changes
- **Terminal motif screening**: Analyze existing contig ends for telomeric sequences
- **Quality metrics**: Comprehensive validation scores for extension confidence

#### Modern Logging System

- **Rich console output**: Beautiful, informative terminal displays using the Rich library
- **Flexible log levels**: Configurable verbosity (DEBUG, INFO, WARNING, ERROR, CRITICAL)
- **File logging support**: Optional detailed log file output with timestamps
- **Context-aware messaging**: Function-specific logging with proper error context

#### Enhanced File I/O

- **Multiple output formats**: Support for both FASTA and FASTQ output with quality scores
- **Flexible organization**: Multi-file output with automatic directory creation
- **Statistics integration**: Optional overhang statistics in FASTA headers
- **Cross-platform paths**: Robust path handling using pathlib

### 🔧 Infrastructure Improvements

#### Comprehensive Test Suite

- **65 new unit tests**: Complete coverage across 4 new test modules
- **Edge case handling**: Robust testing of error conditions and malformed inputs
- **Mock-based isolation**: Proper test isolation using unittest.mock for external dependencies
- **Integration testing**: End-to-end workflow validation

#### Modern Code Architecture

- **Full type annotations**: Complete type hints throughout codebase for better maintainability
- **Dataclass integration**: Modern Python patterns for structured data handling
- **Exception handling**: Comprehensive error handling with user-friendly messages
- **NumPy-style documentation**: Extensive docstrings following established conventions

#### Updated Development Environment

- **Modern dependency management**: Updated to latest stable versions of core dependencies
- **Pre-commit integration**: Automated code quality checks with enhanced hook configurations
- **Ruff formatting**: Modern Python formatter and linter replacing Black
- **Enhanced CI/CD**: Updated GitHub Actions supporting Python 3.8-3.13

### 🔄 Breaking Changes

#### CLI Interface Changes

- **Command structure**: Changed from `teloclip [options] file` to `teloclip {filter|extract|extend} [options] files`
- **Extract behavior**: Now creates separate files per contig end by default for better organization
- **Output formatting**: New Rich-based logging provides different (improved) output formatting

#### Migration Required

```bash
# Old filter usage
teloclip --ref-idx ref.fa.fai input.sam > output.sam

# New filter usage
teloclip filter --ref-idx ref.fa.fai input.sam > output.sam

# Old extract usage
teloclip-extract --ref-idx ref.fa.fai --extract-reads input.sam

# New extract usage
teloclip extract --ref-idx ref.fa.fai input.sam
```

### 📦 Dependencies

#### Added

- **biopython**: Sequence manipulation and file format support
- **click**: Modern CLI framework with enhanced help and validation
- **pyfaidx**: Efficient FASTA file indexing and random access
- **pysam**: BAM/SAM file processing with C-speed performance
- **rich**: Beautiful terminal output and logging enhancements

#### Development Dependencies Added

- **ruff**: Modern Python linting and formatting
- **pre-commit**: Automated code quality enforcement
- Enhanced GitHub Actions workflows

### 🐛 Bug Fixes

- **SAM processing**: Fixed edge cases in soft-clip extraction and alignment validation
- **Memory management**: Resolved memory leaks in large file processing
- **Path handling**: Improved cross-platform file path management
- **Error reporting**: Enhanced error messages with actionable guidance

### 📚 Documentation

- **Comprehensive README**: Complete rewrite with usage examples for all sub-commands
- **Installation guide**: Updated instructions for new dependency requirements
- **Contributing guidelines**: Enhanced development setup with modern tooling
- **API documentation**: Complete docstring coverage for all public functions

### 🚀 Performance

- **Streaming processing**: Constant memory usage regardless of genome size
- **Compiled regex patterns**: Faster motif matching through pre-compilation
- **Indexed file access**: Dramatic speed improvements for large file processing
- **Efficient data structures**: Optimized memory layout using dataclasses and type hints

### 🎯 Quality Assurance

- **100% test pass rate**: All 65 new unit tests passing across Python 3.8-3.13
- **Cross-platform testing**: Verified compatibility on Linux, macOS, and Windows
- **Memory profiling**: Validated constant memory usage on multi-GB test files
- **Integration validation**: End-to-end testing of complete workflows

---

## Migration Guide for v0.3.0

### For Existing Filter Users (Most Common)

The core filtering functionality remains identical, just add the `filter` sub-command:

```bash
# Before v0.3.0
teloclip --ref-idx genome.fa.fai alignments.bam > filtered.sam

# v0.3.0 and later
teloclip filter --ref-idx genome.fa.fai alignments.bam > filtered.sam
```

### For Existing Extract Users

The extraction functionality is enhanced but requires the new sub-command:

```bash
# Before v0.3.0
teloclip-extract --ref-idx genome.fa.fai --extract-dir output/ alignments.bam

# v0.3.0 and later
teloclip extract --ref-idx genome.fa.fai --extract-dir output/ alignments.bam
```

### New Extension Capability

The major new feature for automatic genome completion:

```bash
# New in v0.3.0 - automatic contig extension
teloclip extend overhangs.bam genome.fa --output-fasta extended.fa --stats-report report.txt

# Preview changes without applying
teloclip extend overhangs.bam genome.fa --dry-run --stats-report preview.txt
```

---

## [Previous Releases]

### [0.2.x] - Previous versions

- Basic alignment filtering functionality
- Simple sequence extraction
- Command-line interface for telomere analysis

### [0.1.x] - Initial releases

- Core soft-clip detection algorithms
- Basic motif matching capabilities
- Foundation SAM/BAM processing tools

---

**Note**: This major release (v0.3.0) represents the largest enhancement to Teloclip since its creation, transforming it from a specialized filtering tool into a comprehensive platform for telomere analysis and genome completion. While breaking changes are introduced in the CLI interface, the core functionality remains fully compatible with simple command updates.
