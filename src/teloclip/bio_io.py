"""
Memory-efficient FASTA I/O using BioPython SeqIO.

This module provides streaming FASTA reading and writing functions that don't
load entire genomes into memory, making it suitable for large genome processing.
"""

from pathlib import Path
import sys
from typing import Dict, Iterator, Tuple, Union

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def stream_fasta_sequences(
    fasta_path: Union[str, Path],
) -> Iterator[Tuple[str, str, str]]:
    """
    Stream FASTA sequences without loading entire file into memory.

    Parameters
    ----------
    fasta_path : Union[str, Path]
        Path to FASTA file.

    Yields
    ------
    Tuple[str, str, str]
        (sequence_name, description, sequence) for each sequence in the file.
    """
    with open(fasta_path, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            yield record.id, record.description, str(record.seq)


def load_fasta_sequences(fasta_path: Union[str, Path]) -> Dict[str, Tuple[str, str]]:
    """
    Load FASTA sequences into memory (compatible with existing fasta2dict format).

    This function maintains compatibility with existing code while using BioPython
    for parsing. Returns the same format as fasta2dict: {name: (header, sequence)}.

    Parameters
    ----------
    fasta_path : Union[str, Path]
        Path to FASTA file.

    Returns
    -------
    Dict[str, Tuple[str, str]]
        Dictionary mapping sequence names to (description, sequence) tuples.
    """
    sequences = {}
    with open(fasta_path, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            sequences[record.id] = (record.description, str(record.seq))
    return sequences


def write_fasta_sequences(
    sequences: Dict[str, Tuple[str, str]],
    output_path: Union[str, Path, None] = None,
    line_length: int = 80,
) -> None:
    """
    Write sequences to FASTA file or stdout in a memory-efficient manner.

    Parameters
    ----------
    sequences : Dict[str, Tuple[str, str]]
        Dictionary mapping sequence names to (description, sequence) tuples.
    output_path : Union[str, Path, None], optional
        Output file path. If None, writes to stdout.
    line_length : int, optional
        Number of characters per line in output FASTA (default: 80).
    """
    # Create SeqRecord objects
    seq_records = []
    for seq_id, (description, sequence) in sequences.items():
        record = SeqRecord(Seq(sequence), id=seq_id, description=description)
        seq_records.append(record)

    # Write to file or stdout
    if output_path is None:
        # Write to stdout
        for record in seq_records:
            SeqIO.write(record, sys.stdout, 'fasta')
    else:
        with open(output_path, 'w') as handle:
            for record in seq_records:
                SeqIO.write(record, handle, 'fasta')


def stream_write_fasta_sequences(
    sequences: Iterator[Tuple[str, str, str]],
    output_path: Union[str, Path, None] = None,
) -> None:
    """
    Write sequences from iterator to FASTA file or stdout.

    This function is truly memory-efficient as it processes sequences one at a time
    without loading them all into memory.

    Parameters
    ----------
    sequences : Iterator[Tuple[str, str, str]]
        Iterator yielding (seq_id, description, sequence) tuples.
    output_path : Union[str, Path, None], optional
        Output file path. If None, writes to stdout.
    """

    def create_records():
        """
        Convert sequence tuples to SeqRecord objects for BioPython output.

        Yields
        ------
        SeqRecord
            BioPython SeqRecord object for each input sequence tuple.
        """
        for seq_id, description, sequence in sequences:
            yield SeqRecord(Seq(sequence), id=seq_id, description=description)

    if output_path is None:
        # Write to stdout
        SeqIO.write(create_records(), sys.stdout, 'fasta')
    else:
        with open(output_path, 'w') as handle:
            SeqIO.write(create_records(), handle, 'fasta')


def validate_fasta_against_fai(
    fasta_path: Union[str, Path], fai_dict: Dict[str, int]
) -> Tuple[set, set]:
    """
    Validate that FASTA file contains sequences referenced in FAI index.

    Parameters
    ----------
    fasta_path : Union[str, Path]
        Path to FASTA file.
    fai_dict : Dict[str, int]
        Dictionary from read_fai with sequence names and lengths.

    Returns
    -------
    Tuple[set, set]
        (missing_from_fasta, missing_from_fai) - sets of sequence names that are
        missing from FASTA file or missing from FAI index respectively.
    """
    fasta_sequences = set()

    # Get sequence names from FASTA
    with open(fasta_path, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            fasta_sequences.add(record.id)

    fai_sequences = set(fai_dict.keys())

    missing_from_fasta = fai_sequences - fasta_sequences
    missing_from_fai = fasta_sequences - fai_sequences

    return missing_from_fasta, missing_from_fai
