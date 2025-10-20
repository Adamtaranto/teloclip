"""
Enhanced I/O module for streaming telomere sequence extraction.

This module provides memory-efficient FASTA/FASTQ writing with buffering,
statistics tracking, and BioPython integration for the extract command.
"""

from collections import defaultdict
from pathlib import Path
import sys
from typing import Dict, Optional, Union

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class ExtractionStats:
    """Track extraction statistics and generate reports."""

    def __init__(self):
        """Initialize statistics counters."""
        self.total_sam_lines = 0
        self.filter_counts = {'soft_clip': 0, 'quality': 0, 'anchor': 0}
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
