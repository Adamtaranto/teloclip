"""
I/O module for streaming telomere sequence extraction.

This module provides memory-efficient FASTA/FASTQ writing with buffering,
statistics tracking, and BioPython integration for the extract command.
"""

from collections import defaultdict
import logging
from pathlib import Path
import sys
from typing import Dict, Optional, Union

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class EfficientSequenceWriter:
    """
    Memory-efficient sequence writer with buffering and BioPython integration.

    Parameters
    ----------
    output_path : Union[str, Path, None], optional
        Output file path. If None, writes to stdout. Default is None.
    output_format : str, optional
        Output format ('fasta' or 'fastq'). Default is 'fasta'.
    buffer_size : int, optional
        Number of sequences to buffer before writing. Default is 1000.
    """

    def __init__(
        self,
        output_path: Optional[Union[str, Path]] = None,
        output_format: str = 'fasta',
        buffer_size: int = 1000,
    ):
        """Initialize writer with optional buffering.

        Parameters
        ----------
        output_path : Union[str, Path, None]
            Output file path. If None, writes to stdout.
        output_format : str
            Output format ('fasta' or 'fastq').
        buffer_size : int
            Number of sequences to buffer before writing.
        """
        self.output_path = Path(output_path) if output_path else None
        self.output_format = output_format.lower()
        self.buffer_size = buffer_size
        self.buffer = []
        self.file_handle = None
        self.sequences_written = 0

        # Validate format
        if self.output_format not in {'fasta', 'fastq'}:
            raise ValueError(f'Unsupported output format: {output_format}')

    def __enter__(self):
        """
        Context manager entry.

        Returns
        -------
        EfficientSequenceWriter
            Returns self for context manager usage.
        """
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Context manager exit.

        Parameters
        ----------
        exc_type : type
            Exception type if an exception occurred.
        exc_val : Exception
            Exception value if an exception occurred.
        exc_tb : traceback
            Exception traceback if an exception occurred.
        """
        self.close()

    def open(self):
        """
        Open output file handle.
        """
        if self.output_path is None:
            self.file_handle = sys.stdout
        else:
            # Create parent directories if needed
            self.output_path.parent.mkdir(parents=True, exist_ok=True)
            self.file_handle = open(self.output_path, 'w')

    def write_sequence(
        self,
        seq_id: str,
        sequence: str,
        description: str = '',
        quality: Optional[str] = None,
        stats: Optional[Dict] = None,
    ):
        """
        Write sequence with optional stats in header.

        Parameters
        ----------
        seq_id : str
            Sequence identifier.
        sequence : str
            DNA sequence.
        description : str, optional
            Sequence description.
        quality : str, optional
            Quality string (required for FASTQ format).
        stats : dict, optional
            Statistics to include in header.
        """
        # Build enhanced description with stats
        if stats:
            stat_parts = []
            if 'mapq' in stats:
                stat_parts.append(f'mapq={stats["mapq"]}')
            if 'clip_length' in stats:
                stat_parts.append(f'clip_len={stats["clip_length"]}')
            if 'overhang_length' in stats:
                stat_parts.append(f'overhang_len={stats["overhang_length"]}')
            if 'motif_counts' in stats and stats['motif_counts']:
                motif_str = ','.join(
                    f'{k}:{v}' for k, v in stats['motif_counts'].items()
                )
                stat_parts.append(f'motifs={motif_str}')

            if stat_parts and description:
                description = f'{description} | {" ".join(stat_parts)}'
            elif stat_parts:
                description = ' | '.join(stat_parts)

        # Create SeqRecord
        record = SeqRecord(Seq(sequence), id=seq_id, description=description)

        # Add quality for FASTQ
        if self.output_format == 'fastq':
            if quality is None:
                # Generate dummy quality if not provided
                quality = 'I' * len(sequence)  # Phred+33 quality of 40
            record.letter_annotations['phred_quality'] = [ord(q) - 33 for q in quality]

        # Add to buffer
        self.buffer.append(record)

        # Flush if buffer is full
        if len(self.buffer) >= self.buffer_size:
            self.flush()

    def flush(self):
        """
        Force write buffered sequences.

        Returns
        -------
        None
            No return value.
        """
        if self.buffer and self.file_handle:
            SeqIO.write(self.buffer, self.file_handle, self.output_format)
            self.sequences_written += len(self.buffer)
            self.buffer.clear()

    def close(self):
        """
        Close writer and cleanup.
        """
        # Write any remaining buffered sequences
        if self.buffer:
            self.flush()

        # Close file handle if not stdout
        if self.file_handle and self.file_handle != sys.stdout:
            self.file_handle.close()
            self.file_handle = None


class MultiFileSequenceWriter:
    """
    Context manager for writing sequences to multiple files organized by contig and end type.

    This class creates separate output files for each contig and end type combination
    (e.g., contig1_L.fasta, contig1_R.fasta) and manages file handles efficiently.

    Parameters
    ----------
    base_dir : str or Path
        Base directory where output files will be created.
    prefix : str, optional
        Prefix for output filenames. Default is "teloclip".
    output_format : str, optional
        Output format, either 'fasta' or 'fastq'. Default is 'fasta'.
    buffer_size : int, optional
        Buffer size for file writing. Default is 8192.
    use_sam_attributes : bool, optional
        Format statistics as SAM attributes for FASTQ output. Default is False.

    Examples
    --------
    >>> with MultiFileSequenceWriter('/output', 'sample') as writer:
    ...     writer.write_sequence('chr1', 'L', 'read1', 'ATCG', 'description')
    """

    def __init__(
        self,
        base_dir: Union[str, Path],
        prefix: str = 'teloclip',
        output_format: str = 'fasta',
        buffer_size: int = 8192,
        use_sam_attributes: bool = False,
    ):
        """Initialize the multi-file sequence writer."""
        self.base_dir = Path(base_dir)
        self.prefix = prefix
        self.output_format = output_format.lower()
        self.buffer_size = buffer_size
        self.use_sam_attributes = use_sam_attributes
        self.file_handles = {}
        self.sequence_counts = defaultdict(int)

        # Validate output format
        if self.output_format not in ('fasta', 'fastq'):
            raise ValueError(f'Unsupported output format: {self.output_format}')

    def __enter__(self):
        """
        Enter context manager and create base directory.

        Returns
        -------
        MultiFileSequenceWriter
            Returns self for use in with statement.
        """
        self.base_dir.mkdir(parents=True, exist_ok=True)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Exit context manager and close all file handles.

        Parameters
        ----------
        exc_type : type or None
            Exception type if an exception was raised, None otherwise.
        exc_val : Exception or None
            Exception value if an exception was raised, None otherwise.
        exc_tb : traceback or None
            Exception traceback if an exception was raised, None otherwise.
        """
        for handle in self.file_handles.values():
            handle.close()
        self.file_handles.clear()

    def _get_file_handle(self, contig_name: str, end: str):
        """
        Get or create file handle for a specific contig and end combination.

        Parameters
        ----------
        contig_name : str
            Name of the contig.
        end : str
            End type ('L' for left, 'R' for right).

        Returns
        -------
        file handle
            Open file handle for writing.
        """
        key = (contig_name, end)

        if key not in self.file_handles:
            # Create filename
            if self.prefix:
                filename = f'{self.prefix}_{contig_name}_{end}.{self.output_format}'
            else:
                filename = f'{contig_name}_{end}.{self.output_format}'
            filepath = self.base_dir / filename

            # Open file handle
            self.file_handles[key] = open(filepath, 'w', buffering=self.buffer_size)

        return self.file_handles[key]

    def write_sequence(
        self,
        contig_name: str,
        end: str,
        seq_id: str,
        sequence: str,
        description: str = '',
        stats: Optional[Dict] = None,
    ):
        """
        Write a sequence to the appropriate file.

        Parameters
        ----------
        contig_name : str
            Name of the contig.
        end : str
            End type ('L' for left, 'R' for right).
        seq_id : str
            Sequence identifier.
        sequence : str
            DNA sequence.
        description : str, optional
            Sequence description. Default is empty string.
        stats : dict, optional
            Statistics to include in sequence header.

        Raises
        ------
        ValueError
            If output format is not supported or required parameters are missing.
        """
        if not contig_name or not end or not seq_id or not sequence:
            raise ValueError('contig_name, end, seq_id, and sequence are required')

        # Get file handle
        handle = self._get_file_handle(contig_name, end)

        # Build header
        header = seq_id
        if description:
            header += f' {description}'

        # Add stats to header if provided
        if stats:
            if self.use_sam_attributes and self.output_format == 'fastq':
                # Format as SAM attributes: tag:type:value
                stats_parts = []
                for k, v in stats.items():
                    if isinstance(v, int):
                        stats_parts.append(f'{k}:i:{v}')
                    elif isinstance(v, float):
                        stats_parts.append(f'{k}:f:{v}')
                    else:
                        stats_parts.append(f'{k}:Z:{v}')
                stats_str = ' '.join(stats_parts)
            else:
                # Format as key=value pairs (default for FASTA or when SAM attributes disabled)
                stats_str = ' '.join(f'{k}={v}' for k, v in stats.items())
            header += f' {stats_str}'

        if self.output_format == 'fasta':
            # Write FASTA format
            handle.write(f'>{header}\n')
            # Write sequence in lines of 80 characters
            for i in range(0, len(sequence), 80):
                handle.write(f'{sequence[i : i + 80]}\n')

        elif self.output_format == 'fastq':
            # Write FASTQ format (with dummy quality scores)
            quality = 'I' * len(sequence)  # Phred+33 quality score of 40
            handle.write(f'@{header}\n')
            handle.write(f'{sequence}\n')
            handle.write('+\n')
            handle.write(f'{quality}\n')

        # Update counts
        key = (contig_name, end)
        self.sequence_counts[key] += 1

    def get_sequence_counts(self) -> Dict[tuple, int]:
        """
        Get the count of sequences written for each contig-end combination.

        Returns
        -------
        dict
            Dictionary mapping (contig_name, end) tuples to sequence counts.
        """
        return dict(self.sequence_counts)


class ExtractionStats:
    """Track extraction statistics and generate reports."""

    def __init__(self):
        """Initialize statistics counters."""
        self.total_sam_lines = 0
        self.filter_counts = {
            'unmapped': 0,
            'secondary': 0,
            'soft_clip': 0,
            'quality': 0,
            'anchor': 0,
            'max_break': 0,
            'min_clip': 0,
            'motifs': 0,
        }
        self.total_filtered = 0
        self.total_passed_to_split = 0
        # Track which contigs have overhangs found
        self.contigs_with_overhangs = set()

        # Legacy fields for backward compatibility
        self.total_alignments = 0
        self.left_overhangs = 0
        self.right_overhangs = 0
        self.contigs_processed = set()  # Contigs with valid overhangs
        self.contigs_with_left = set()
        self.contigs_with_right = set()
        self.motif_matches = defaultdict(int)
        self.quality_filtered = 0
        self.anchor_filtered = 0
        self.soft_clip_filtered = 0  # Lines without soft clips
        self.malformed_lines = 0

    def record_sam_line(self):
        """Record that we processed a SAM line."""
        self.total_sam_lines += 1

    def record_filter(self, filter_type):
        """
        Record that an alignment was filtered out.

        Parameters
        ----------
        filter_type : str
            Type of filter applied ('soft_clip', 'quality', or 'anchor').
        """
        if filter_type in self.filter_counts:
            self.filter_counts[filter_type] += 1
            self.total_filtered += 1

    def record_passed_to_split(self):
        """Record that an alignment passed filtering and went to split processing."""
        self.total_passed_to_split += 1

    def record_contig_with_overhang(self, contig_name):
        """
        Record that this contig had at least one overhang processed.

        Parameters
        ----------
        contig_name : str
            Name of the contig that had overhangs found.
        """
        self.contigs_with_overhangs.add(contig_name)

    def record_alignment(
        self,
        contig_name: str,
        is_left: bool,
        motif_counts: Optional[Dict[str, int]] = None,
    ):
        """
        Record processed alignment statistics.

        Parameters
        ----------
        contig_name : str
            Name of the contig being processed.
        is_left : bool
            Whether this is a left-end overhang alignment.
        motif_counts : dict, optional
            Dictionary of motif patterns and their counts.
        """
        self.total_alignments += 1
        self.contigs_processed.add(contig_name)
        # Also record for the new tracking
        self.contigs_with_overhangs.add(contig_name)

        if is_left:
            self.left_overhangs += 1
            self.contigs_with_left.add(contig_name)
        else:
            self.right_overhangs += 1
            self.contigs_with_right.add(contig_name)

        if motif_counts:
            for motif, count in motif_counts.items():
                self.motif_matches[motif] += count

    def generate_report(self, reference_contigs: Optional[set] = None) -> str:
        """
        Generate a comprehensive extraction report.

        Parameters
        ----------
        reference_contigs : set, optional
            Set of all contig names from reference index.

        Returns
        -------
        str
            Formatted report string.
        """
        lines = []
        lines.append('Extraction Statistics Report')
        lines.append('=' * 40)

        # SAM processing summary
        lines.append(f'Total SAM lines processed: {self.total_sam_lines:,}')
        lines.append(
            f'Lines passed to split processing: {self.total_passed_to_split:,}'
        )
        lines.append(f'Total lines filtered: {self.total_filtered:,}')

        # Filtering breakdown
        if self.total_filtered > 0:
            lines.append('\nFiltering breakdown:')
            for filter_type, count in self.filter_counts.items():
                if count > 0:
                    percentage = (
                        (count / self.total_sam_lines) * 100
                        if self.total_sam_lines > 0
                        else 0
                    )
                    lines.append(f'  {filter_type}: {count:,} ({percentage:.1f}%)')

        # Contig analysis
        if reference_contigs:
            total_ref_contigs = len(reference_contigs)
            contigs_with_overhangs = len(self.contigs_with_overhangs)
            contigs_with_zero_overhangs = total_ref_contigs - contigs_with_overhangs

            lines.append('\nContig Analysis:')
            lines.append(f'Total reference contigs: {total_ref_contigs:,}')
            lines.append(f'Contigs with overhangs: {contigs_with_overhangs:,}')
            lines.append(
                f'Contigs with zero overhangs: {contigs_with_zero_overhangs:,}'
            )

            if contigs_with_zero_overhangs > 0:
                zero_overhang_contigs = reference_contigs - self.contigs_with_overhangs
                lines.append(
                    f'Zero-overhang contigs: {", ".join(sorted(zero_overhang_contigs))}'
                )

        # Legacy alignment statistics
        lines.append(f'\nTotal alignments processed: {self.total_alignments:,}')
        lines.append(f'Left overhangs: {self.left_overhangs:,}')
        lines.append(f'Right overhangs: {self.right_overhangs:,}')

        return '\n'.join(lines)

    def log_exclusion_summary(self):
        """Log comprehensive exclusion criteria summary."""
        logging.info(
            f'Processed {self.total_sam_lines} SAM records.\n'
            f'Passed to processing: {self.total_passed_to_split}\n'
            f'Exclusion summary:\n'
            f'  - Unmapped reads: {self.filter_counts["unmapped"]}\n'
            f'  - Secondary alignments: {self.filter_counts["secondary"]}\n'
            f'  - No soft-clips: {self.filter_counts["soft_clip"]}\n'
            f'  - Below quality threshold: {self.filter_counts["quality"]}\n'
            f'  - Below min_anchor threshold: {self.filter_counts["anchor"]}\n'
            f'  - Beyond max_break threshold: {self.filter_counts["max_break"]}\n'
            f'  - Below min_clip threshold: {self.filter_counts["min_clip"]}\n'
            f'  - No telomeric motifs: {self.filter_counts["motifs"]}\n'
            f'Total filtered: {self.total_filtered} alignments after all filtering.'
        )


class BufferedContigWriter:
    """
    Memory-efficient writer for contig sequences with buffering.

    Provides buffered writing to FASTA/FASTQ files with automatic flushing
    when buffer size limits are reached.

    Parameters
    ----------
    output_path : str or Path
        Output file path.
    buffer_size : int, default 1000
        Number of sequences to buffer before writing.
    file_format : str, default "fasta"
        Output format ("fasta" or "fastq").
    """

    def __init__(
        self,
        output_path: Union[str, Path],
        buffer_size: int = 1000,
        file_format: str = 'fasta',
    ):
        """
        Initialize the buffered writer.

        Parameters
        ----------
        output_path : str or Path
            Output file path.
        buffer_size : int, default 1000
            Number of sequences to buffer before writing.
        file_format : str, default "fasta"
            Output format ("fasta" or "fastq").
        """
        self.output_path = Path(output_path)
        self.buffer_size = buffer_size
        self.file_format = file_format.lower()
        self.buffer = []
        self.total_written = 0

        if self.file_format not in ['fasta', 'fastq']:
            raise ValueError("file_format must be 'fasta' or 'fastq'")

    def add_sequence(
        self,
        sequence_id: str,
        sequence: str,
        quality: Optional[str] = None,
        description: str = '',
    ):
        """
        Add a sequence to the buffer.

        Parameters
        ----------
        sequence_id : str
            Sequence identifier.
        sequence : str
            Sequence string.
        quality : str, optional
            Quality scores for FASTQ format.
        description : str, default ""
            Sequence description.
        """
        if self.file_format == 'fastq' and quality is None:
            raise ValueError('Quality scores required for FASTQ format')

        record = SeqRecord(Seq(sequence), id=sequence_id, description=description)

        if self.file_format == 'fastq' and quality:
            record.letter_annotations['phred_quality'] = [ord(c) - 33 for c in quality]

        self.buffer.append(record)

        if len(self.buffer) >= self.buffer_size:
            self.flush()

    def flush(self):
        """
        Write buffered sequences to file.

        Returns
        -------
        None
            Returns nothing.
        """
        if not self.buffer:
            return

        mode = 'a' if self.total_written > 0 else 'w'

        with open(self.output_path, mode) as handle:
            SeqIO.write(self.buffer, handle, self.file_format)

        self.total_written += len(self.buffer)
        self.buffer.clear()

    def close(self):
        """
        Flush remaining sequences and close writer.

        Returns
        -------
        None
            Returns nothing.
        """
        self.flush()

    def __enter__(self):
        """
        Context manager entry.

        Returns
        -------
        BufferedContigWriter
            Returns self for use in with statement.
        """
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Context manager exit.

        Parameters
        ----------
        exc_type : type
            Exception type.
        exc_val : Exception
            Exception value.
        exc_tb : traceback
            Exception traceback.
        """
        self.close()


class StreamingGenomeProcessor:
    """
    Memory-efficient genome processing with indexed access.

    Uses pysam for indexed FASTA access to avoid loading entire genomes
    into memory while maintaining random access capabilities.

    Parameters
    ----------
    fasta_path : str or Path
        Path to indexed FASTA file (.fai index required).
    """

    def __init__(self, fasta_path: Union[str, Path]):
        """
        Initialize the streaming processor.

        Parameters
        ----------
        fasta_path : str or Path
            Path to indexed FASTA file (.fai index required).
        """
        try:
            import pysam
        except ImportError as exc:
            raise ImportError(
                'pysam is required for streaming genome processing'
            ) from exc

        self.fasta_path = Path(fasta_path)
        self.fasta = pysam.FastaFile(str(self.fasta_path))

        # Verify index exists
        index_path = Path(str(self.fasta_path) + '.fai')
        if not index_path.exists():
            raise FileNotFoundError(
                f'FASTA index not found: {index_path}. '
                "Create with 'samtools faidx' or pysam.faidx()"
            )

    def get_sequence(
        self, contig: str, start: Optional[int] = None, end: Optional[int] = None
    ) -> str:
        """
        Extract sequence from specified genomic region.

        Parameters
        ----------
        contig : str
            Contig/chromosome name.
        start : int, optional
            Start position (0-based, inclusive).
        end : int, optional
            End position (0-based, exclusive).

        Returns
        -------
        str
            Extracted sequence.
        """
        try:
            if start is None and end is None:
                return self.fasta.fetch(contig)
            elif start is None:
                return self.fasta.fetch(contig, 0, end)
            elif end is None:
                return self.fasta.fetch(contig, start)
            else:
                return self.fasta.fetch(contig, start, end)
        except KeyError as exc:
            raise KeyError(f"Contig '{contig}' not found in FASTA file") from exc
        except Exception as e:
            raise RuntimeError(f'Error fetching sequence: {e}') from e

    def get_contig_length(self, contig: str) -> int:
        """
        Get the length of a contig.

        Parameters
        ----------
        contig : str
            Contig name.

        Returns
        -------
        int
            Contig length in bases.
        """
        try:
            return self.fasta.get_reference_length(contig)
        except KeyError as exc:
            raise KeyError(f"Contig '{contig}' not found in FASTA file") from exc

    def list_contigs(self) -> list:
        """
        Get list of all contig names in the FASTA file.

        Returns
        -------
        list
            List of contig names.
        """
        return list(self.fasta.references)

    def close(self):
        """
        Close the FASTA file handle.

        Returns
        -------
        None
            Returns nothing.
        """
        if hasattr(self, 'fasta'):
            self.fasta.close()

    def __enter__(self):
        """
        Context manager entry.

        Returns
        -------
        StreamingGenomeProcessor
            Returns self for use in with statement.
        """
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Context manager exit.

        Parameters
        ----------
        exc_type : type
            Exception type.
        exc_val : Exception
            Exception value.
        exc_tb : traceback
            Exception traceback.
        """
        self.close()


def write_sequences_streaming(
    sequences: Dict[str, Dict[str, str]],
    output_path: Union[str, Path],
    file_format: str = 'fasta',
    buffer_size: int = 1000,
) -> int:
    """
    Write sequences to file using streaming I/O.

    Parameters
    ----------
    sequences : dict
        Nested dictionary: {contig: {side: sequence}}.
    output_path : str or Path
        Output file path.
    file_format : str, default "fasta"
        Output format ("fasta" or "fastq").
    buffer_size : int, default 1000
        Buffer size for writing.

    Returns
    -------
    int
        Number of sequences written.
    """
    total_written = 0

    with BufferedContigWriter(
        output_path, buffer_size=buffer_size, file_format=file_format
    ) as writer:
        for contig_name, sides in sequences.items():
            for side, sequence in sides.items():
                if sequence:  # Skip empty sequences
                    seq_id = f'{contig_name}_{side}'
                    writer.add_sequence(seq_id, sequence)
                    total_written += 1

    return total_written


def create_fasta_index(fasta_path: Union[str, Path]) -> Path:
    """
    Create FASTA index file if it doesn't exist.

    Parameters
    ----------
    fasta_path : str or Path
        Path to FASTA file.

    Returns
    -------
    Path
        Path to created index file.
    """
    try:
        import pysam
    except ImportError as exc:
        raise ImportError('pysam is required for FASTA indexing') from exc

    fasta_path = Path(fasta_path)
    index_path = Path(str(fasta_path) + '.fai')

    if not index_path.exists():
        print(f'Creating FASTA index: {index_path}', file=sys.stderr)
        pysam.faidx(str(fasta_path))

    return index_path
