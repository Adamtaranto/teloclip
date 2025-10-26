#!/usr/bin/env python3
"""
Generate synthetic test datasets for teloclip extend integration testing.

This script creates synthetic genome contigs and corresponding SAM alignments
with specific characteristics to test all edge cases in the teloclip extend functionality.
"""

import argparse
from pathlib import Path
import random
from typing import Dict


def generate_random_sequence(length: int, gc_content: float = 0.5) -> str:
    """
    Generate a random DNA sequence with specified GC content.

    Parameters
    ----------
    length : int
        Length of sequence to generate.
    gc_content : float
        GC content (0.0 to 1.0). Default is 0.5.

    Returns
    -------
    str
        Random DNA sequence.
    """
    bases = ['A', 'T', 'G', 'C']

    # Calculate base frequencies
    gc_bases = int(length * gc_content)
    at_bases = length - gc_bases

    # Create base pool
    sequence_bases = (
        ['G'] * (gc_bases // 2)
        + ['C'] * (gc_bases // 2)
        + ['A'] * (at_bases // 2)
        + ['T'] * (at_bases // 2)
    )

    # Add remaining bases if odd numbers
    while len(sequence_bases) < length:
        sequence_bases.append(random.choice(bases))

    # Shuffle to randomize order
    random.shuffle(sequence_bases)

    return ''.join(sequence_bases)


def generate_sequence_with_motifs(
    length: int,
    left_motif: str = '',
    right_motif: str = '',
    extension_region: str = '',
    gc_content: float = 0.5,
) -> str:
    """
    Generate DNA sequence with specific motifs at ends and extension regions.

    Parameters
    ----------
    length : int
        Total length of main sequence.
    left_motif : str
        Motif to include near left end.
    right_motif : str
        Motif to include near right end.
    extension_region : str
        Special sequence for overhang extension testing.
    gc_content : float
        GC content for random regions.

    Returns
    -------
    str
        Generated DNA sequence with motifs.
    """
    # Reserve space for motifs
    motif_buffer = 50  # Buffer space around motifs
    core_length = length - len(left_motif) - len(right_motif) - (2 * motif_buffer)

    if core_length < 100:
        core_length = 100
        length = core_length + len(left_motif) + len(right_motif) + (2 * motif_buffer)

    # Generate core sequence
    core_seq = generate_random_sequence(core_length, gc_content)

    # Build sequence with motifs
    sequence_parts = []

    # Left buffer + motif
    if left_motif:
        left_buffer = generate_random_sequence(motif_buffer, gc_content)
        sequence_parts.extend([left_buffer, left_motif])
    else:
        sequence_parts.append(generate_random_sequence(motif_buffer, gc_content))

    # Core sequence
    sequence_parts.append(core_seq)

    # Right motif + buffer
    if right_motif:
        right_buffer = generate_random_sequence(motif_buffer, gc_content)
        sequence_parts.extend([right_motif, right_buffer])
    else:
        sequence_parts.append(generate_random_sequence(motif_buffer, gc_content))

    sequence = ''.join(sequence_parts)

    # Add extension region for testing (this will be in overhang areas)
    if extension_region:
        # This will be added to reads as overhang, not to the contig itself
        pass

    return sequence


def generate_test_contigs() -> Dict[str, str]:
    """
    Generate synthetic contigs with specific characteristics for testing.

    Returns
    -------
    Dict[str, str]
        Dictionary mapping contig names to sequences.
    """
    # Set random seed for reproducibility
    random.seed(42)

    contigs = {}

    # contig_1: Normal contig with moderate overhangs (telomeric ends)
    contigs['contig_1'] = generate_sequence_with_motifs(
        length=1000,
        left_motif='TTAGGGTTAGGG',  # 2x TTAGGG
        right_motif='CCCTAACCCTAA',  # 2x CCCTAA
        gc_content=0.4,
    )

    # contig_2: High-overhang contig (10x more alignments target)
    contigs['contig_2'] = generate_sequence_with_motifs(
        length=2000,
        left_motif='TTAGGGTTAGGGTTAGGG',  # 3x TTAGGG
        right_motif='CCCTAACCCTAACCCTAA',  # 3x CCCTAA
        gc_content=0.45,
    )

    # contig_3: One-sided high overhangs (left end only)
    contigs['contig_3'] = generate_sequence_with_motifs(
        length=1500,
        left_motif='TTAGGGTTAGGGTTAGGGTTAGGG',  # 4x TTAGGG (high signal)
        right_motif='',  # No special motif on right
        gc_content=0.5,
    )

    # contig_4: Homopolymer test contig
    # Main sequence is normal, homopolymer will be in overhang extensions
    contigs['contig_4'] = generate_sequence_with_motifs(
        length=800,
        left_motif='ATGCATGC',  # Simple motif
        right_motif='GCATGCAT',  # Simple motif
        gc_content=0.5,
    )

    # contig_5: Motif-rich contig (TTAGGG repeats)
    contigs['contig_5'] = generate_sequence_with_motifs(
        length=1200,
        left_motif='TTAGGG' * 5,  # 5x TTAGGG
        right_motif='CCCTAA' * 4,  # 4x CCCTAA
        gc_content=0.3,  # AT-rich like telomeres
    )

    # contig_6: No overhangs (control) - completely normal sequence
    contigs['contig_6'] = generate_random_sequence(900, gc_content=0.5)

    # contig_7: Break alignment test contig
    contigs['contig_7'] = generate_sequence_with_motifs(
        length=1100,
        left_motif='GGAATTCC',  # Unique motif for alignment break testing
        right_motif='CCTTAAGG',  # Unique motif for alignment break testing
        gc_content=0.5,
    )

    return contigs


def write_fasta_file(contigs: Dict[str, str], output_path: Path) -> None:
    """
    Write contigs to FASTA file.

    Parameters
    ----------
    contigs : Dict[str, str]
        Dictionary mapping contig names to sequences.
    output_path : Path
        Output FASTA file path.
    """
    with open(output_path, 'w') as f:
        for name, sequence in contigs.items():
            f.write(f'>{name}\n')
            # Write sequence in 80-character lines
            for i in range(0, len(sequence), 80):
                f.write(f'{sequence[i : i + 80]}\n')

    print(f'Generated {len(contigs)} contigs in {output_path}')

    # Print contig summary
    print('\nContig Summary:')
    print('Name\t\tLength\tDescription')
    print('-' * 50)
    descriptions = {
        'contig_1': 'Normal with telomeric ends',
        'contig_2': 'High-coverage target (10x)',
        'contig_3': 'One-sided overhangs (left)',
        'contig_4': 'Homopolymer test target',
        'contig_5': 'Motif-rich (TTAGGG/CCCTAA)',
        'contig_6': 'Control (no special features)',
        'contig_7': 'Break alignment test',
    }

    for name, sequence in contigs.items():
        desc = descriptions.get(name, 'Unknown')
        print(f'{name}\t{len(sequence)}\t{desc}')


def main():
    """Main function to generate test contigs."""
    parser = argparse.ArgumentParser(
        description='Generate synthetic contigs for teloclip testing'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('./test_data'),
        help='Output directory for test files (default: ./test_data)',
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed for reproducibility (default: 42)',
    )

    args = parser.parse_args()

    # Set random seed
    random.seed(args.seed)

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Generate contigs
    print('Generating synthetic test contigs...')
    contigs = generate_test_contigs()

    # Write FASTA file
    fasta_path = args.output_dir / 'synthetic_contigs.fasta'
    write_fasta_file(contigs, fasta_path)

    print('\nTest contigs generated successfully!')
    print(f'FASTA file: {fasta_path}')
    print('\nNext steps:')
    print('1. Index the FASTA file: samtools faidx synthetic_contigs.fasta')
    print('2. Run generate_sam_alignments.py to create corresponding alignments')


if __name__ == '__main__':
    main()
