import pytest
from src.teloclip.motifs import (
    count_continuous_runs,
    construct_regex_pattern,
    format_pattern_counts,
    count_patterns_in_sequence,
    check_sequence_for_patterns,
)


def test_check_sequence_for_patterns():
    # Test case with a simple DNA sequence and patterns
    dna_sequence = "ATCGATCGATCG"
    regex_patterns = ["ATC", "GAT", "CGA"]
    assert check_sequence_for_patterns(dna_sequence, regex_patterns) is True

    # Test case with an empty DNA sequence
    assert check_sequence_for_patterns("", ["ATC", "GAT"]) is False

    # Test case with an empty pattern list
    assert check_sequence_for_patterns("ATCG", []) is False

    # Test case with no matches in the DNA sequence
    dna_sequence = "ATCGATCGATCG"
    regex_patterns = ["AAA", "GGG", "CCC"]
    assert check_sequence_for_patterns(dna_sequence, regex_patterns) is False

    # Test case with a complex DNA sequence and patterns
    dna_sequence = "ATCGGATCGAGCGCGAATCG"
    regex_patterns = ["ATC", "GAT", "CGA"]
    assert check_sequence_for_patterns(dna_sequence, regex_patterns) is True


def test_count_patterns_in_sequence():
    # Test case with a simple DNA sequence and patterns
    dna_sequence = "ATCGATCGATCG"
    regex_patterns = ["ATC", "GAT", "CGA"]
    assert count_patterns_in_sequence(dna_sequence, regex_patterns) == {
        "ATC": 3,
        "GAT": 2,
        "CGA": 2,
    }

    # Test case with an empty DNA sequence
    assert count_patterns_in_sequence("", ["ATC", "GAT"]) == {"ATC": 0, "GAT": 0}

    # Test case with an empty pattern list
    assert count_patterns_in_sequence("ATCG", []) == {}

    # Test case with a complex DNA sequence and patterns
    dna_sequence = "ATCGGATCGAGCGCGAATCG"
    regex_patterns = ["ATC", "GAT", "CGA"]
    assert count_patterns_in_sequence(dna_sequence, regex_patterns) == {
        "ATC": 3,
        "GAT": 1,
        "CGA": 2,
    }


def test_format_pattern_counts():
    # Test case with a simple pattern counts dictionary
    pattern_counts = {"ATC": 3, "GAT": 1, "CGA": 2}
    assert format_pattern_counts(pattern_counts) == "ATC=3:GAT=1:CGA=2"

    # Test case with an empty pattern counts dictionary
    assert format_pattern_counts({}) == ""

    # Test case with a single pattern count
    assert format_pattern_counts({"ATC": 1}) == "ATC=1"

    # Test case with patterns and counts containing special characters
    pattern_counts = {"A(TC)": 2, "G+AT": 3, "C.A": 1}
    assert format_pattern_counts(pattern_counts) == "A(TC)=2:G+AT=3:C.A=1"


def test_count_continuous_runs():
    # Test case with a simple sequence
    assert count_continuous_runs("AAATTTGGCCCG") == [
        ("A", 3),
        ("T", 3),
        ("G", 2),
        ("C", 3),
        ("G", 1),
    ]

    # Test case with an empty string
    assert count_continuous_runs("") == []

    # Test case with a single character
    assert count_continuous_runs("A") == [("A", 1)]

    # Test case with alternating characters
    assert count_continuous_runs("ATATAT") == [
        ("A", 1),
        ("T", 1),
        ("A", 1),
        ("T", 1),
        ("A", 1),
        ("T", 1),
    ]


def test_construct_regex_pattern():
    # Test case with a simple motif
    motif_sequence = "TTTAGGG"
    motif_tuples = count_continuous_runs(motif_sequence)
    assert construct_regex_pattern(motif_tuples) == r"T{2,4}AG{2,4}"

    # Test case with a motif of length 1
    motif_sequence = "A"
    motif_tuples = count_continuous_runs(motif_sequence)
    assert construct_regex_pattern(motif_tuples) == r"A"

    # Test case with an empty motif
    assert construct_regex_pattern([]) == ""

    # Test case with a motif of alternating characters
    motif_sequence = "ATATAT"
    motif_tuples = count_continuous_runs(motif_sequence)
    assert construct_regex_pattern(motif_tuples) == r"ATATAT"

    # Additional test case with a more complex motif
    motif_sequence = "AAAAGGGTTTCCCC"
    motif_tuples = count_continuous_runs(motif_sequence)
    assert construct_regex_pattern(motif_tuples) == r"A{3,5}G{2,4}T{2,4}C{3,5}"


# Run the tests
if __name__ == "__main__":
    pytest.main()
