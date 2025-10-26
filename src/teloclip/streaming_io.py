"""
Memory-efficient streaming I/O for large genome processing.

This module provides streaming functions for processing large genomes without
loading entire files into memory, using pysam for indexed access.
"""

from pathlib import Path
import sys
from typing import Dict, Iterator, Optional, Tuple, Union

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam


class StreamingGenomeProcessor:
    """
    Memory-efficient processor for large genomes using indexed access.

    This class provides methods to process genomes contig-by-contig without
    loading the entire genome into memory.

    Parameters
    ----------
    fasta_path : Union[str, Path]
        Path to indexed FASTA file (.fai index must exist).
    bam_path : Union[str, Path]
        Path to indexed BAM file (.bai index must exist).
    """

    def __init__(self, fasta_path: Union[str, Path], bam_path: Union[str, Path]):
        """
        Initialize the processor with indexed files.

        Parameters
        ----------
        fasta_path : Union[str, Path]
            Path to indexed FASTA file (.fai index must exist).
        bam_path : Union[str, Path]
            Path to indexed BAM file (.bai index must exist).
        """
        self.fasta_path = Path(fasta_path)
        self.bam_path = Path(bam_path)

        # Verify indexes exist
        fai_path = Path(str(self.fasta_path) + '.fai')
        bai_path = Path(str(self.bam_path) + '.bai')

        if not fai_path.exists():
            raise FileNotFoundError(f'FASTA index not found: {fai_path}')
        if not bai_path.exists():
            raise FileNotFoundError(f'BAM index not found: {bai_path}')

        # Initialize file handles (will be opened when needed)
        self._fasta_file = None
        self._bam_file = None

    def __enter__(self):
        """
        Context manager entry.

        Returns
        -------
        StreamingGenomeProcessor
            Self instance for context management.
        """
        self._fasta_file = pysam.FastaFile(str(self.fasta_path))
        self._bam_file = pysam.AlignmentFile(str(self.bam_path), 'rb')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Context manager exit.

        Parameters
        ----------
        exc_type : type or None
            Exception type if an exception occurred.
        exc_val : Exception or None
            Exception instance if an exception occurred.
        exc_tb : traceback or None
            Traceback object if an exception occurred.
        """
        if self._fasta_file:
            self._fasta_file.close()
        if self._bam_file:
            self._bam_file.close()

    def get_contig_sequence(self, contig_name: str) -> str:
        """
        Load sequence for a single contig.

        Parameters
        ----------
        contig_name : str
            Name of the contig to load.

        Returns
        -------
        str
            The contig sequence.
        """
        if not self._fasta_file:
            raise RuntimeError('Processor not initialized. Use as context manager.')

        try:
            return self._fasta_file.fetch(contig_name)
        except KeyError as exc:
            raise KeyError(f"Contig '{contig_name}' not found in FASTA file") from exc

    def get_contig_alignments(self, contig_name: str) -> Iterator:
        """
        Get alignments for a specific contig.

        Parameters
        ----------
        contig_name : str
            Name of the contig.

        Yields
        ------
        pysam.AlignedSegment
            Alignment records for the specified contig.
        """
        if not self._bam_file:
            raise RuntimeError('Processor not initialized. Use as context manager.')

        try:
            for alignment in self._bam_file.fetch(contig_name):
                yield alignment
        except ValueError:
            # Contig not found in BAM - return empty iterator
            return

    def get_all_contig_names(self) -> Iterator[str]:
        """
        Get names of all contigs in the FASTA file.

        Yields
        ------
        str
            Contig names.
        """
        if not self._fasta_file:
            raise RuntimeError('Processor not initialized. Use as context manager.')

        for contig_name in self._fasta_file.references:
            yield contig_name


class BufferedContigWriter:
    """
    Buffered writer for extended contigs to minimize I/O operations.

    Parameters
    ----------
    output_path : Optional[Union[str, Path]]
        Output file path. If None, writes to stdout.
    buffer_size : int
        Buffer size in bytes (default: 1MB).
    """

    def __init__(
        self,
        output_path: Optional[Union[str, Path]] = None,
        buffer_size: int = 1024 * 1024,
    ):
        """
        Initialize buffered writer.

        Parameters
        ----------
        output_path : Optional[Union[str, Path]]
            Output file path. If None, writes to stdout.
        buffer_size : int
            Buffer size in bytes (default: 1MB).
        """
        self.output_path = Path(output_path) if output_path else None
        self.buffer_size = buffer_size
        self.buffer = []
        self.current_size = 0
        self._file_handle = None

    def __enter__(self):
        """
        Context manager entry.

        Returns
        -------
        BufferedContigWriter
            Self instance for context management.
        """
        if self.output_path:
            self._file_handle = open(self.output_path, 'w')
        else:
            self._file_handle = sys.stdout
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Context manager exit.

        Parameters
        ----------
        exc_type : type or None
            Exception type if an exception occurred.
        exc_val : Exception or None
            Exception instance if an exception occurred.
        exc_tb : traceback or None
            Traceback object if an exception occurred.
        """
        self.flush()
        if self.output_path and self._file_handle:
            self._file_handle.close()

    def write_contig(self, name: str, sequence: str, description: str = ''):
        """
        Write a contig to the buffer.

        Parameters
        ----------
        name : str
            Contig name.
        sequence : str
            Contig sequence.
        description : str, optional
            Contig description.
        """
        # Create SeqRecord and format as FASTA
        record = SeqRecord(Seq(sequence), id=name, description=description)
        fasta_str = record.format('fasta')

        self.buffer.append(fasta_str)
        self.current_size += len(fasta_str)

        # Flush if buffer is full
        if self.current_size >= self.buffer_size:
            self.flush()

    def flush(self):
        """Flush buffer to output."""
        if self.buffer and self._file_handle:
            for fasta_record in self.buffer:
                self._file_handle.write(fasta_record)
            self._file_handle.flush()

            # Clear buffer
            self.buffer.clear()
            self.current_size = 0


def stream_sam_lines_for_contig(
    bam_file: pysam.AlignmentFile, contig_name: str
) -> Iterator[str]:
    """
    Stream SAM format lines for a specific contig.

    Parameters
    ----------
    bam_file : pysam.AlignmentFile
        Opened BAM file.
    contig_name : str
        Name of the contig to process.

    Yields
    ------
    str
        SAM format lines for alignments to the specified contig.
    """
    try:
        for alignment in bam_file.fetch(contig_name):
            # Skip unmapped reads
            if alignment.is_unmapped:
                continue
            # Skip secondary/supplementary alignments
            if alignment.is_secondary or alignment.is_supplementary:
                continue

            # Convert pysam alignment back to SAM line format
            sam_line = alignment.to_string()
            yield sam_line
    except ValueError:
        # Contig not found in BAM
        return


def validate_indexed_files(
    fasta_path: Union[str, Path], bam_path: Union[str, Path]
) -> Tuple[bool, str]:
    """
    Validate that required index files exist.

    Parameters
    ----------
    fasta_path : Union[str, Path]
        Path to FASTA file.
    bam_path : Union[str, Path]
        Path to BAM file.

    Returns
    -------
    Tuple[bool, str]
        (is_valid, error_message).
    """
    fasta_path = Path(fasta_path)
    bam_path = Path(bam_path)

    # Check FASTA file exists
    if not fasta_path.exists():
        return False, f'FASTA file not found: {fasta_path}'

    # Check BAM file exists
    if not bam_path.exists():
        return False, f'BAM file not found: {bam_path}'

    # Check FASTA index exists
    fai_path = Path(str(fasta_path) + '.fai')
    if not fai_path.exists():
        return (
            False,
            f'FASTA index not found: {fai_path}. Create with: samtools faidx {fasta_path}',
        )

    # Check BAM index exists
    bai_path = Path(str(bam_path) + '.bai')
    if not bai_path.exists():
        return (
            False,
            f'BAM index not found: {bai_path}. Create with: samtools index {bam_path}',
        )

    return True, ''


def copy_unmodified_contigs(
    processor: StreamingGenomeProcessor,
    writer: BufferedContigWriter,
    modified_contigs: set,
    contig_dict: Dict[str, int],
) -> None:
    """
    Copy unmodified contigs from input to output.

    Parameters
    ----------
    processor : StreamingGenomeProcessor
        Processor for accessing original sequences.
    writer : BufferedContigWriter
        Writer for output.
    modified_contigs : set
        Set of contig names that were modified (to skip).
    contig_dict : Dict[str, int]
        Dictionary of all contigs and their lengths.
    """
    for contig_name in contig_dict:
        if contig_name not in modified_contigs:
            try:
                sequence = processor.get_contig_sequence(contig_name)
                writer.write_contig(contig_name, sequence)
            except KeyError:
                # Contig not found in FASTA (should have been caught earlier)
                continue
