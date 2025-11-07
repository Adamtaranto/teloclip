"""
Motif analysis functions for sequence pattern matching.

This module provides functions for exact and fuzzy motif matching in DNA sequences,
including regex pattern generation, motif counting, and analysis utilities.
"""

import re
from typing import Dict, List, Tuple


def make_motif_regex(motif: str) -> str:
    """
    Create a regex pattern to match exact motifs.

    Parameters
    ----------
    motif : str
        The input motif string to create a regex pattern for.

    Returns
    -------
    str
        The regex pattern for exact motif matching.
    """
    # Escape special characters in the motif
    escaped_motif = re.escape(motif)
    # Return the regex pattern for the exact motif
    return rf'({escaped_motif})'


def make_fuzzy_motif_regex(motif: str) -> str:
    """
    Create a regex pattern to match fuzzy motifs with runs of characters
    that differ by plus or minus one compared to the original motif.

    Parameters
    ----------
    motif : str
        The input motif string.

    Returns
    -------
    str
        The constructed regex pattern as a raw string.
    """
    # Count continuous runs of characters in the motif
    motif_tuples = count_continuous_runs(motif)
    # Construct the regex pattern based on the counted runs
    pattern = construct_regex_pattern(motif_tuples)
    # Return the final regex pattern allowing for the specified minimum repeats
    return rf'({pattern})'


def count_regex_patterns_in_sequence(
    dna_sequence: str, regex_patterns: List[str]
) -> Dict[str, int]:
    """
    Count occurrences of regular expression patterns in a DNA sequence.

    Parameters
    ----------
    dna_sequence : str
        The input DNA sequence.
    regex_patterns : List[str]
        List of regular expression patterns to count.

    Returns
    -------
    Dict[str, int]
        Dictionary mapping pattern names to occurrence counts.
    """
    # Confirm that regex patterns does not end with min repeats quantifier of format {x,y} or {x,}
    for pattern in regex_patterns:
        if re.search(r'\{\d+,\d+\}$', pattern) or re.search(r'\{\d+,\}$', pattern):
            raise ValueError(
                f"Pattern '{pattern}' ends with min repeats quantifier. "
                'Please provide patterns that do not end with min repeats quantifiers.'
            )

    pattern_counts = {}  # Initialize an empty dictionary to store counts

    # Iterate through the list of regex patterns
    for pattern in regex_patterns:
        # Use re.findall to find all occurrences of the pattern in the DNA sequence
        matches = re.findall(pattern, dna_sequence)
        # Count the number of matches and store in the dictionary
        pattern_counts[pattern] = len(matches)

    # Or as a dict comprehension
    # pattern_counts = {pattern: len(re.findall(pattern, dna_sequence)) for pattern in regex_patterns}

    return pattern_counts


def count_continuous_runs(dna_string: str) -> list:
    """
    Count the length of all continuous runs of a DNA base in a string.

    Parameters
    ----------
    dna_string : str
        The input DNA string.

    Returns
    -------
    list
        A list of tuples where each tuple contains a character and the number
        of times it occurred consecutively in the input DNA string.
    """

    # Check if the input string is empty
    if not dna_string:
        return []

    # Initialize variables
    current_base = dna_string[0]
    current_count = 1
    result = []

    # Iterate through the DNA string starting from the second character
    for base in dna_string[1:]:
        if base == current_base:
            # If the current base is the same as the previous one, increment the count
            current_count += 1
        else:
            # If the current base is different, append the tuple and reset the count
            result.append((current_base, current_count))
            current_base = base
            current_count = 1

    # Append the last tuple after the loop
    result.append((current_base, current_count))

    return result


def construct_regex_pattern(motif_tuples: List[Tuple[str, int]]) -> str:
    """
    Construct a regex pattern to match sequences with runs of characters
    that differ by plus or minus one compared to the original input sequence.

    Parameters
    ----------
    motif_tuples : List[Tuple[str, int]]
        List of tuples where each tuple contains a character and the number
        of times it occurred consecutively.

    Returns
    -------
    str
        The constructed regex pattern as a raw string.
    """

    pattern_parts = []

    # Note: This idea was adapted from https://github.com/JanaSperschneider/FindTelomeres/blob/master/FindTelomeres.py

    for char, count in motif_tuples:
        if count == 1:
            pattern_parts.append(
                re.escape(char)
            )  # If count is 1, just escape the character
        else:
            # If count is greater than 1, add a range allowing for plus or minus one
            pattern_parts.append(rf'{re.escape(char)}{{{count - 1},{count + 1}}}')

    return rf'{"".join(pattern_parts)}'


def format_pattern_counts(pattern_counts: Dict[str, int]) -> str:
    """
    Format a dictionary of pattern counts into a string.

    Parameters
    ----------
    pattern_counts : Dict[str, int]
        A dictionary where keys are regex patterns and values are counts
        of each pattern.

    Returns
    -------
    str
        A formatted string representation of the pattern counts.
    """
    # Use a list comprehension to generate formatted pairs of pattern=count
    formatted_pairs = [
        f'{pattern}={count}' for pattern, count in pattern_counts.items()
    ]

    # Join the formatted pairs with ":" and return the result
    return ':'.join(formatted_pairs)


def check_sequence_for_patterns(
    dna_sequence: str, regex_patterns: List[str], min_repeats: int = 1
) -> bool:
    """
    Check a DNA sequence for instances of one or more regular expressions.

    Parameters
    ----------
    dna_sequence : str
        The input DNA sequence.
    regex_patterns : List[str]
        List of regular expression patterns to check against.
    min_repeats : int, optional
        Minimum number of sequential repeats required for a match.
        Default is 1.

    Returns
    -------
    bool
        True if any of the patterns match the sequence, False otherwise.
    """
    # Modify patterns to include minimum repeats requirement
    if min_repeats > 1:
        regex_patterns = [
            rf'({pattern}){{{min_repeats},}}' for pattern in regex_patterns
        ]
    # Use any() to check if any pattern matches the DNA sequence
    return any(re.search(pattern, dna_sequence) for pattern in regex_patterns)
