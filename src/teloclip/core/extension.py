"""
Contig extension algorithms and utilities.

This module provides functionality to extend draft contigs using selected
overhang sequences from soft-clipped alignments.
"""

from typing import Tuple

from .analysis import OverhangInfo


def calculate_extension_position(
    alignment_pos: int, alignment_end: int, contig_length: int, is_left: bool
) -> Tuple[int, int]:
    """
    Calculate the precise position for extending a contig and required trimming.

    Parameters
    ----------
    alignment_pos : int
        Start position of alignment (1-based).
    alignment_end : int
        End position of alignment (1-based, inclusive).
    contig_length : int
        Length of the contig.
    is_left : bool
        True if this is a left overhang, False for right.

    Returns
    -------
    Tuple[int, int]
        (extension_position, trim_length) where:
        - extension_position: where to add the overhang (0-based)
        - trim_length: how much to trim from contig end (0 if no trimming needed)
    """
    if is_left:
        # For left extensions, we extend at position 0
        # Trim if alignment doesn't start exactly at position 1
        trim_length = max(0, alignment_pos - 1)
        return 0, trim_length
    else:
        # For right extensions, we extend at the end
        # Trim if alignment doesn't end exactly at contig end
        trim_length = max(0, contig_length - alignment_end)
        extension_position = contig_length - trim_length
        return extension_position, trim_length


def trim_contig_end(sequence: str, trim_length: int, is_left_end: bool) -> str:
    """
    Trim bases from the end of a contig sequence.

    Parameters
    ----------
    sequence : str
        Original contig sequence.
    trim_length : int
        Number of bases to trim.
    is_left_end : bool
        True to trim from left end, False to trim from right end.

    Returns
    -------
    str
        Trimmed sequence.
    """
    if trim_length <= 0:
        return sequence

    if is_left_end:
        return sequence[trim_length:]
    else:
        return sequence[:-trim_length]


def extend_contig(
    sequence: str, overhang: OverhangInfo, position: int, is_left_end: bool
) -> str:
    """
    Extend a contig sequence with an overhang sequence.

    Parameters
    ----------
    sequence : str
        Original contig sequence.
    overhang : OverhangInfo
        Overhang information containing the sequence to add.
    position : int
        Position to insert the overhang (0-based).
    is_left_end : bool
        True if extending left end, False if extending right end.

    Returns
    -------
    str
        Extended sequence.
    """
    if is_left_end:
        # Add overhang to the beginning
        return overhang.sequence + sequence
    else:
        # Add overhang to the end
        return sequence + overhang.sequence


def validate_extension(
    original: str, extended: str, overhang: OverhangInfo, trim_length: int = 0
) -> bool:
    """
    Validate that a contig extension was performed correctly.

    The length comparison is made against the *untrimmed* original sequence,
    which is what makes a net-shortening extension detectable. The containment
    check is necessarily made against the trimmed remainder, since the trimmed
    bases are no longer present in the extended sequence.

    Parameters
    ----------
    original : str
        Original, untrimmed contig sequence.
    extended : str
        Extended contig sequence.
    overhang : OverhangInfo
        Overhang that was added.
    trim_length : int, optional
        Number of bases trimmed from the extended end before the overhang was
        grafted on (default: 0).

    Returns
    -------
    bool
        True if extension is valid, False otherwise.
    """
    # The part of the original sequence that survives trimming, and so must
    # still be present in the extended sequence.
    if trim_length > 0:
        retained = (
            original[trim_length:] if overhang.is_left else original[:-trim_length]
        )
    else:
        retained = original

    if overhang.is_left:
        # Check that overhang was added to the beginning
        if not extended.startswith(overhang.sequence):
            return False
        # Check that the retained portion of the contig survived
        if retained not in extended:
            return False
    else:
        # Check that overhang was added to the end
        if not extended.endswith(overhang.sequence):
            return False
        # Check that the retained portion of the contig survived
        if retained not in extended:
            return False

    # The extension must leave the contig longer than it started. Comparing
    # against the untrimmed original is essential: comparing against the
    # trimmed remainder makes this check trivially true and lets an extension
    # that removes more bases than it adds pass silently.
    return len(extended) > len(original)


def apply_contig_extension(
    contig_seq: str,
    overhang: OverhangInfo,
    contig_length: int,
    min_net_gain: int = 1,
) -> Tuple[str, dict]:
    """
    Apply a complete contig extension including position calculation and trimming.

    Extending a contig end trims the bases the supporting read did not cover and
    grafts the read's whole soft clip on in their place, so the net change in
    length is the number of clipped bases lying past the contig terminus. This
    function refuses to apply an extension whose net change falls below
    ``min_net_gain``, which is what stops a short clip anchored well inside the
    contig from silently truncating it.

    Parameters
    ----------
    contig_seq : str
        Original contig sequence.
    overhang : OverhangInfo
        Overhang to use for extension.
    contig_length : int
        Length of the original contig.
    min_net_gain : int, optional
        Minimum acceptable increase in contig length, in bases (default: 1).

    Returns
    -------
    Tuple[str, dict]
        (extended_sequence, extension_info) where extension_info contains:
        - 'overhang_length': length of added overhang
        - 'net_gain': net change in contig length (overhang minus trim)
        - 'trim_length': number of bases trimmed
        - 'extension_position': where extension was applied
        - 'original_length': original contig length
        - 'final_length': final contig length
        - 'read_name': source read name

    Raises
    ------
    ValueError
        If the extension would grow the contig by fewer than ``min_net_gain``
        bases, or if the resulting sequence fails validation.
    """
    # Calculate extension position and trimming
    position, trim_length = calculate_extension_position(
        overhang.alignment_pos, overhang.alignment_end, contig_length, overhang.is_left
    )

    # Apply trimming if needed
    working_seq = trim_contig_end(contig_seq, trim_length, overhang.is_left)

    # Apply extension
    extended_seq = extend_contig(working_seq, overhang, position, overhang.is_left)

    # Reject any extension that does not grow the contig. Every byte of
    # sequence mutation flows through this function, so this is the one place
    # the guarantee can be enforced for both ends and both call paths.
    net_gain = len(extended_seq) - len(contig_seq)
    if net_gain < min_net_gain:
        side = 'left' if overhang.is_left else 'right'
        raise ValueError(
            f'Extension from read {overhang.read_name} on the {side} end of '
            f'contig {overhang.contig_name} would change its length by '
            f'{net_gain:+d} bp (trim {trim_length}, add {overhang.length}); '
            f'a net gain of at least {min_net_gain} bp is required'
        )

    # Validate extension against the untrimmed original
    if not validate_extension(contig_seq, extended_seq, overhang, trim_length):
        raise ValueError(f'Extension validation failed for read {overhang.read_name}')

    # Prepare extension information
    extension_info = {
        'overhang_length': overhang.length,
        'net_gain': net_gain,
        'trim_length': trim_length,
        'extension_position': position,
        'original_length': len(contig_seq),
        'final_length': len(extended_seq),
        'read_name': overhang.read_name,
        'is_left': overhang.is_left,
    }

    return extended_seq, extension_info
