import re
from typing import List, Tuple, Dict


def format_pattern_counts(pattern_counts: Dict[str, int]) -> str:
    """
    Format a dictionary of pattern counts into a string.

    Parameters:
    - pattern_counts (Dict[str, int]): A dictionary where keys are regex patterns
                                    and values are counts of each pattern.

    Returns:
    - str: A formatted string with pattern counts.
    """
    # Use a list comprehension to generate formatted pairs of pattern=count
    formatted_pairs = [
        f"{pattern}={count}" for pattern, count in pattern_counts.items()
    ]

    # Join the formatted pairs with ":" and return the result
    return ":".join(formatted_pairs)


def count_patterns_in_sequence(
    dna_sequence: str, regex_patterns: List[str]
) -> Dict[str, int]:
    """
    Count occurrences of regular expression patterns in a DNA sequence.

    Parameters:
    - dna_sequence (str): The input DNA sequence.
    - regex_patterns (List[str]): List of regular expression patterns to count.

    Returns:
    - Dict[str, int]: A dictionary where keys are regex patterns and values are
                     counts of each pattern in the sequence.
    """
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


def check_sequence_for_patterns(dna_sequence: str, regex_patterns: List[str]) -> bool:
    """
    Check a DNA sequence for instances of one or more regular expressions.

    Parameters:
    - dna_sequence (str): The input DNA sequence.
    - regex_patterns (List[str]): List of regular expression patterns to check against.

    Returns:
    - bool: True if any of the patterns match the sequence, False otherwise.
    """
    # Use any() to check if any pattern matches the DNA sequence
    return any(re.search(pattern, dna_sequence) for pattern in regex_patterns)


def count_continuous_runs(dna_string: str) -> list:
    """
    Count the length of all continuous runs of a DNA base in a string.

    Parameters:
    - dna_string (str): The input DNA string.

    Returns:
    - list: A list of tuples where each tuple contains a character and the number
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

    Parameters:
    - motif_tuples (List[Tuple[str, int]]): List of tuples where each tuple contains
                                            a character and the number of times it
                                            occurred consecutively.

    Returns:
    - str: The constructed regex pattern as a raw string.
    """

    pattern_parts = []

    for char, count in motif_tuples:
        if count == 1:
            pattern_parts.append(
                re.escape(char)
            )  # If count is 1, just escape the character
        else:
            # If count is greater than 1, add a range allowing for plus or minus one
            pattern_parts.append(rf"{re.escape(char)}{{{count-1},{count+1}}}")

    return rf"{''.join(pattern_parts)}"
