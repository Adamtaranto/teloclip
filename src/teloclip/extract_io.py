"""
Enhanced I/O module for efficient sequence extraction.

This module provides memory-efficient FASTA/FASTQ writing with buffering,
statistics tracking, and BioPython integration for the extract command.
"""

import sys
from pathlib import Path
from typing import Dict, Optional, Union
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class ExtractionStats:
    """Track extraction statistics and generate reports."""

    def __init__(self):
        """Initialize statistics counters."""
        self.total_alignments = 0
        self.left_overhangs = 0
        self.right_overhangs = 0
        self.contigs_processed = set()
        self.contigs_with_left = set()
        self.contigs_with_right = set()
        self.motif_matches = defaultdict(int)
        self.quality_filtered = 0
        self.anchor_filtered = 0

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

        if is_left:
            self.left_overhangs += 1
            self.contigs_with_left.add(contig_name)
        else:
            self.right_overhangs += 1
            self.contigs_with_right.add(contig_name)

        if motif_counts:
            for motif, count in motif_counts.items():
                self.motif_matches[motif] += count

    def record_filter(self, filter_type: str):
        """
        Record filtered alignments.

        Parameters
        ----------
        filter_type : str
            Type of filter applied ('quality' or 'anchor').
        """
        if filter_type == 'quality':
            self.quality_filtered += 1
        elif filter_type == 'anchor':
            self.anchor_filtered += 1

    def generate_report(self) -> str:
        """
        Generate comprehensive extraction report.

        Returns
        -------
        str
            Formatted markdown report with extraction statistics.
        """

        lines = []
        lines.append('# Teloclip Extract Statistics Report')
        lines.append('=' * 50)
        lines.append('')

        lines.append('## Processing Summary')
        lines.append(f'Total alignments processed: {self.total_alignments}')
        lines.append(f'Left end overhangs: {self.left_overhangs}')
        lines.append(f'Right end overhangs: {self.right_overhangs}')
        lines.append('')

        lines.append('## Contig Summary')
        lines.append(f'Total contigs with overhangs: {len(self.contigs_processed)}')
        lines.append(f'Contigs with left overhangs: {len(self.contigs_with_left)}')
        lines.append(f'Contigs with right overhangs: {len(self.contigs_with_right)}')
        lines.append('')

        if self.quality_filtered or self.anchor_filtered:
            lines.append('## Filtering Summary')
            if self.quality_filtered:
                lines.append(f'Quality filtered: {self.quality_filtered}')
            if self.anchor_filtered:
                lines.append(f'Anchor length filtered: {self.anchor_filtered}')
            lines.append('')

        if self.motif_matches:
            lines.append('## Motif Analysis')
            for motif, count in sorted(self.motif_matches.items()):
                lines.append(f'{motif}: {count} matches')
            lines.append('')

        return '\\n'.join(lines)


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
    Manage multiple sequence writers for different contigs/ends.

    Parameters
    ----------
    base_dir : Union[str, Path, None], optional
        Base output directory. Default is None (uses current directory).
    prefix : str, optional
        Prefix for output filenames. Default is None.
    output_format : str, optional
        Output format ('fasta' or 'fastq'). Default is 'fasta'.
    buffer_size : int, optional
        Buffer size for each writer. Default is 1000.
    """

    def __init__(
        self,
        base_dir: Optional[Union[str, Path]] = None,
        prefix: Optional[str] = None,
        output_format: str = 'fasta',
        buffer_size: int = 1000,
    ):
        """Initialize multi-file writer manager.

        Parameters
        ----------
        base_dir : Union[str, Path, None]
            Base output directory.
        prefix : str, optional
            Prefix for output filenames.
        output_format : str
            Output format ('fasta' or 'fastq').
        buffer_size : int
            Buffer size for each writer.
        """
        self.base_dir = Path(base_dir) if base_dir else Path.cwd()
        self.prefix = prefix
        self.output_format = output_format
        self.buffer_size = buffer_size
        self.writers: Dict[str, EfficientSequenceWriter] = {}

        # Create output directory
        self.base_dir.mkdir(parents=True, exist_ok=True)

    def get_writer(self, contig_name: str, end: str) -> EfficientSequenceWriter:
        """
        Get or create writer for specific contig and end.

        Parameters
        ----------
        contig_name : str
            Name of contig.
        end : str
            End type ('L' for left, 'R' for right).

        Returns
        -------
        EfficientSequenceWriter
            Writer for this contig/end combination.
        """

        key = f'{contig_name}_{end}'

        if key not in self.writers:
            # Build filename
            if self.prefix:
                filename = f'{self.prefix}_{contig_name}_{end}.{self.output_format}'
            else:
                filename = f'{contig_name}_{end}.{self.output_format}'

            output_path = self.base_dir / filename

            # Create and open writer
            writer = EfficientSequenceWriter(
                output_path=output_path,
                output_format=self.output_format,
                buffer_size=self.buffer_size,
            )
            writer.open()
            self.writers[key] = writer

        return self.writers[key]

    def write_sequence(
        self,
        contig_name: str,
        end: str,
        seq_id: str,
        sequence: str,
        description: str = '',
        quality: Optional[str] = None,
        stats: Optional[Dict] = None,
    ):
        """
        Write sequence to appropriate file.

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
            Sequence description.
        quality : str, optional
            Quality string for FASTQ format.
        stats : dict, optional
            Statistics to include in header.
        """
        writer = self.get_writer(contig_name, end)
        writer.write_sequence(seq_id, sequence, description, quality, stats)

    def close_all(self):
        """
        Close all writers.

        Returns
        -------
        None
            No return value.
        """
        for writer in self.writers.values():
            writer.close()
        self.writers.clear()

    def __enter__(self):
        """
        Context manager entry.

        Returns
        -------
        MultiFileSequenceWriter
            Returns self for context manager usage.
        """
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
        self.close_all()
