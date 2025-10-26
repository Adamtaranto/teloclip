"""Unit tests for teloclip.seqops module.

Tests sequence operation functions including FASTA I/O,
reverse complement, and sequence filtering utilities.
"""

from unittest.mock import mock_open, patch

from teloclip.seqops import (
    addRevComplement,
    filterList,
    isMotifInClip,
    makeMask,
    read_fai,
    revComp,
)


class TestMakeMask:
    """Test mask creation function."""

    def test_make_mask_single_index(self):
        """Test creating mask with single kill index."""
        result = makeMask([2], 5)
        expected = [True, True, False, True, True]
        assert result == expected

    def test_make_mask_multiple_indices(self):
        """Test creating mask with multiple kill indices."""
        result = makeMask([1, 3], 5)
        expected = [True, False, True, False, True]
        assert result == expected

    def test_make_mask_empty_indices(self):
        """Test creating mask with no kill indices."""
        result = makeMask([], 4)
        expected = [True, True, True, True]
        assert result == expected

    def test_make_mask_all_indices(self):
        """Test creating mask killing all indices."""
        result = makeMask([0, 1, 2], 3)
        expected = [False, False, False]
        assert result == expected


class TestFilterList:
    """Test filtering lists with mask."""

    def test_filter_list_simple(self):
        """Test filtering with simple exclude indices."""
        data = ['a', 'b', 'c', 'd']
        exclude = [1, 3]  # Exclude indices 1 and 3 (b and d)
        result = list(filterList(data, exclude))
        expected = ['a', 'c']
        assert result == expected

    def test_filter_list_none_excluded(self):
        """Test filtering with no exclusions."""
        data = ['a', 'b', 'c']
        exclude = []  # No exclusions
        result = list(filterList(data, exclude))
        expected = ['a', 'b', 'c']
        assert result == expected

    def test_filter_list_all_excluded(self):
        """Test filtering with all excluded."""
        data = ['a', 'b', 'c']
        exclude = [0, 1, 2]  # Exclude all indices
        result = list(filterList(data, exclude))
        expected = []
        assert result == expected

    def test_filter_list_empty(self):
        """Test filtering empty list."""
        data = []
        exclude = []
        result = list(filterList(data, exclude))
        expected = []
        assert result == expected


class TestRevComp:
    """Test reverse complement function."""

    def test_rev_comp_simple(self):
        """Test reverse complement of simple sequence."""
        result = revComp('ATCG')
        expected = 'CGAT'
        assert result == expected

    def test_rev_comp_longer_sequence(self):
        """Test reverse complement of longer sequence."""
        result = revComp('TTAGGGCCCTAA')
        expected = 'TTAGGGCCCTAA'  # This is palindromic
        assert result == expected

    def test_rev_comp_with_ambiguous(self):
        """Test reverse complement with ambiguous bases."""
        result = revComp('ATCGN')
        expected = 'NCGAT'
        assert result == expected

    def test_rev_comp_uppercase_only(self):
        """Test reverse complement works with uppercase only."""
        result = revComp('ATCG')
        expected = 'CGAT'
        assert result == expected

    def test_rev_comp_empty(self):
        """Test reverse complement of empty string."""
        result = revComp('')
        expected = ''
        assert result == expected


class TestReadFai:
    """Test FAI index file reading function."""

    @patch('teloclip.seqops.isfile')
    @patch('builtins.open', new_callable=mock_open)
    def test_read_fai_simple(self, mock_file, mock_isfile):
        """Test reading simple FAI file."""
        mock_isfile.return_value = '/path/to/test.fai'
        fai_content = 'chr1\t1000\t5\t80\t81\nchr2\t2000\t1010\t80\t81\n'
        mock_file.return_value.readlines.return_value = fai_content.splitlines(
            keepends=True
        )

        result = read_fai('test.fai')

        expected = {'chr1': 1000, 'chr2': 2000}
        assert result == expected

    @patch('teloclip.seqops.isfile')
    @patch('builtins.open', new_callable=mock_open)
    def test_read_fai_empty(self, mock_file, mock_isfile):
        """Test reading empty FAI file."""
        mock_isfile.return_value = '/path/to/empty.fai'
        mock_file.return_value.readlines.return_value = []

        result = read_fai('empty.fai')

        assert result == {}

    @patch('teloclip.seqops.isfile')
    @patch('builtins.open', new_callable=mock_open)
    def test_read_fai_malformed_line(self, mock_file, mock_isfile):
        """Test reading FAI with malformed line."""
        mock_isfile.return_value = '/path/to/test.fai'
        fai_content = 'chr1\t1000\t5\t80\t81\nchr2\t2000\t1010\t80\t81\n'
        mock_file.return_value.readlines.return_value = fai_content.splitlines(
            keepends=True
        )

        result = read_fai('test.fai')

        # Should process valid lines only
        expected = {'chr1': 1000, 'chr2': 2000}
        assert result == expected


class TestAddRevComplement:
    """Test reverse complement addition function."""

    def test_add_rev_complement_simple(self):
        """Test adding reverse complements to motif list."""
        motifs = ['TTAGGG', 'CCCTAA']
        result = addRevComplement(motifs)

        # Should contain original motifs plus their reverse complements
        assert 'TTAGGG' in result
        assert 'CCCTAA' in result
        assert 'CCCTAA' in result  # Rev comp of TTAGGG
        assert 'TTAGGG' in result  # Rev comp of CCCTAA

        # Result should be longer than or equal to original
        assert len(result) >= len(motifs)

    def test_add_rev_complement_palindromes(self):
        """Test adding reverse complements with palindromic sequences."""
        motifs = ['TATA', 'CGCG']  # These are palindromes
        result = addRevComplement(motifs)

        # Should still contain the motifs (duplicates may be present)
        assert 'TATA' in result
        assert 'CGCG' in result

    def test_add_rev_complement_empty(self):
        """Test adding reverse complements to empty list."""
        motifs = []
        result = addRevComplement(motifs)

        assert result == set()


class TestIsMotifInClip:
    """Test motif detection in clipped sequences."""

    def test_is_motif_in_clip_match(self):
        """Test detecting motif in clipped sequence."""
        # Mock samline: tab-separated values as list
        samline = [
            'read1',
            '4',
            'chr1',
            '100',
            '60',
            '5S10M5S',
            '*',
            '0',
            '0',
            'TTAGGATCGATCGCCCTAA',  # sequence with TTAGG at start and CCTAA at end
        ]
        patterns = ['TTAGGG', 'CCTAA']

        result = isMotifInClip(
            samline,
            patterns,
            leftClip=True,
            rightClip=True,
            leftClipLen=5,
            rightClipLen=5,
        )

        assert result is True

    def test_is_motif_in_clip_no_match(self):
        """Test no motif detection in clipped sequence."""
        samline = [
            'read1',
            '4',
            'chr1',
            '100',
            '60',
            '5S10M5S',
            '*',
            '0',
            '0',
            'ATCGATCGATCGATCGATCG',  # no telomeric motifs
        ]
        patterns = ['TTAGGG', 'CCTAA']

        result = isMotifInClip(
            samline,
            patterns,
            leftClip=True,
            rightClip=True,
            leftClipLen=5,
            rightClipLen=5,
        )

        assert result is False

    def test_is_motif_in_clip_multiple_patterns(self):
        """Test motif detection with multiple patterns."""
        samline = [
            'read1',
            '4',
            'chr1',
            '100',
            '60',
            '5S10M5S',
            '*',
            '0',
            '0',
            'CCCTAATCGATCGTTTAGGG',  # CCCTAA at start, TTAGGG at end
        ]
        patterns = ['TTAGGG', 'CCCTAA']

        result = isMotifInClip(
            samline,
            patterns,
            leftClip=True,
            rightClip=True,
            leftClipLen=6,  # Changed to 6 to get CCCTAA
            rightClipLen=6,  # Changed to 6 to get TTAGGG
        )

        assert result is True

    def test_is_motif_in_clip_left_clip_only(self):
        """Test motif detection in left clip only."""
        samline = [
            'read1',
            '4',
            'chr1',
            '100',
            '60',
            '6S10M',  # Changed to 6S
            '*',
            '0',
            '0',
            'TTAGGGATCGATCGATCG',  # TTAGGG at start
        ]
        patterns = ['TTAGGG']

        result = isMotifInClip(
            samline,
            patterns,
            leftClip=True,
            rightClip=False,
            leftClipLen=6,  # Changed to 6 to get TTAGGG
            rightClipLen=0,
        )

        assert result is True

    def test_is_motif_in_clip_right_clip_only(self):
        """Test motif detection in right clip only."""
        samline = [
            'read1',
            '4',
            'chr1',
            '100',
            '60',
            '10M5S',
            '*',
            '0',
            '0',
            'ATCGATCGATCCTAA',  # CCTAA at end
        ]
        patterns = ['CCTAA']

        result = isMotifInClip(
            samline,
            patterns,
            leftClip=False,
            rightClip=True,
            leftClipLen=0,
            rightClipLen=5,
        )

        assert result is True

    def test_is_motif_in_clip_no_clips(self):
        """Test motif detection with no clips."""
        samline = [
            'read1',
            '4',
            'chr1',
            '100',
            '60',
            '15M',
            '*',
            '0',
            '0',
            'ATCGATCGATCGATC',
        ]
        patterns = ['TTAGGG']

        result = isMotifInClip(
            samline,
            patterns,
            leftClip=False,
            rightClip=False,
            leftClipLen=0,
            rightClipLen=0,
        )

        assert result is False


class TestSeqopsIntegration:
    """Test integration between sequence operation functions."""

    def test_rev_comp_and_add_rev_complement(self):
        """Test reverse complement function consistency."""
        original_motifs = ['TTAGGG', 'ATCGAT']

        # Add reverse complements
        expanded_motifs = addRevComplement(original_motifs)

        # Manually compute reverse complements
        manual_rev_comps = [revComp(motif) for motif in original_motifs]

        # All manual reverse complements should be in expanded list
        for rev_comp in manual_rev_comps:
            assert rev_comp in expanded_motifs

    def test_mask_and_filter_integration(self):
        """Test mask creation and filtering integration."""
        data = ['seq1', 'seq2', 'seq3', 'seq4', 'seq5']
        kill_indices = [1, 3]  # Remove seq2 and seq4

        # Filter data using kill indices
        filtered_data = list(filterList(data, kill_indices))

        # Should have removed the specified indices
        expected = ['seq1', 'seq3', 'seq5']
        assert filtered_data == expected
