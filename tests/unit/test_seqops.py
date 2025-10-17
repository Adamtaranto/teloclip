"""Unit tests for teloclip.seqops module.

Tests sequence operation functions including FASTA I/O,
reverse complement, and sequence filtering utilities.
"""

from unittest.mock import mock_open, patch
from teloclip.seqops import (
    makeMask,
    filterList,
    revComp,
    writeClip,
    fasta2dict,
    writefasta,
    read_fai,
    addRevComplement,
    isMotifInClip,
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
    """Test list filtering function."""

    def test_filter_list_simple(self):
        """Test filtering with simple mask."""
        data = ['a', 'b', 'c', 'd']
        exclude = [False, True, False, True]
        result = filterList(data, exclude)
        expected = ['a', 'c']
        assert result == expected

    def test_filter_list_none_excluded(self):
        """Test filtering with no exclusions."""
        data = ['a', 'b', 'c']
        exclude = [False, False, False]
        result = filterList(data, exclude)
        expected = ['a', 'b', 'c']
        assert result == expected

    def test_filter_list_all_excluded(self):
        """Test filtering with all excluded."""
        data = ['a', 'b', 'c']
        exclude = [True, True, True]
        result = filterList(data, exclude)
        expected = []
        assert result == expected

    def test_filter_list_empty(self):
        """Test filtering empty list."""
        data = []
        exclude = []
        result = filterList(data, exclude)
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


class TestWriteClip:
    """Test clip writing function."""

    def test_write_clip_basic(self):
        """Test basic clip writing."""
        result = writeClip(1, 4, 10, 'ATCG', 20)
        # Result should be formatted string with index, sequence, and coordinates
        assert isinstance(result, str)
        assert '0001' in result  # Zero-padded index
        assert 'ATCG' in result  # Sequence

    def test_write_clip_different_padding(self):
        """Test clip writing with different zero padding."""
        result = writeClip(99, 3, 5, 'TTAGGG', 100)
        assert isinstance(result, str)
        assert '099' in result  # Zero-padded to 3 digits
        assert 'TTAGGG' in result


class TestFasta2dict:
    """Test FASTA file reading function."""

    @patch('builtins.open', new_callable=mock_open)
    def test_fasta2dict_simple(self, mock_file):
        """Test reading simple FASTA file."""
        fasta_content = '>seq1\nATCGATCG\n>seq2\nTTAGGGCC\n'
        mock_file.return_value.__iter__.return_value = fasta_content.splitlines()

        result = fasta2dict('test.fasta')

        expected = {'seq1': 'ATCGATCG', 'seq2': 'TTAGGGCC'}
        assert result == expected

    @patch('builtins.open', new_callable=mock_open)
    def test_fasta2dict_multiline_sequences(self, mock_file):
        """Test reading FASTA with multi-line sequences."""
        fasta_content = '>seq1\nATCG\nATCG\n>seq2\nTTAG\nGGCC\n'
        mock_file.return_value.__iter__.return_value = fasta_content.splitlines()

        result = fasta2dict('test.fasta')

        expected = {'seq1': 'ATCGATCG', 'seq2': 'TTAGGGCC'}
        assert result == expected

    @patch('builtins.open', new_callable=mock_open)
    def test_fasta2dict_empty_file(self, mock_file):
        """Test reading empty FASTA file."""
        mock_file.return_value.__iter__.return_value = []

        result = fasta2dict('empty.fasta')

        assert result == {}


class TestWriteFasta:
    """Test FASTA writing function."""

    @patch('builtins.open', new_callable=mock_open)
    def test_write_fasta_simple(self, mock_file):
        """Test writing simple FASTA entry."""
        writefasta('output.fasta', 'test_seq', 'ATCGATCG')

        # Check that file was opened for writing
        mock_file.assert_called_once_with('output.fasta', 'w')

        # Check that proper FASTA format was written
        handle = mock_file.return_value
        written_calls = handle.write.call_args_list

        # Should write header and sequence
        written_text = ''.join(call[0][0] for call in written_calls)
        assert '>test_seq' in written_text
        assert 'ATCGATCG' in written_text

    @patch('builtins.open', new_callable=mock_open)
    def test_write_fasta_long_sequence(self, mock_file):
        """Test writing FASTA with line wrapping."""
        long_seq = 'A' * 100  # 100 A's
        writefasta('output.fasta', 'long_seq', long_seq, length=50)

        handle = mock_file.return_value
        written_calls = handle.write.call_args_list
        written_text = ''.join(call[0][0] for call in written_calls)

        # Should wrap at specified length
        lines = written_text.strip().split('\n')
        seq_lines = [line for line in lines if not line.startswith('>')]

        # All sequence lines except potentially the last should be exactly 50 chars
        for line in seq_lines[:-1]:
            assert len(line) == 50

    @patch('builtins.open', new_callable=mock_open)
    def test_write_fasta_append_mode(self, mock_file):
        """Test writing FASTA in append mode."""
        writefasta('output.fasta', 'test_seq', 'ATCG', append=True)

        # Should open in append mode
        mock_file.assert_called_once_with('output.fasta', 'a')


class TestReadFai:
    """Test FAI index file reading function."""

    @patch('builtins.open', new_callable=mock_open)
    def test_read_fai_simple(self, mock_file):
        """Test reading simple FAI file."""
        fai_content = 'chr1\t1000\t5\t80\t81\nchr2\t2000\t1010\t80\t81\n'
        mock_file.return_value.__iter__.return_value = fai_content.splitlines()

        result = read_fai('test.fai')

        expected = {'chr1': 1000, 'chr2': 2000}
        assert result == expected

    @patch('builtins.open', new_callable=mock_open)
    def test_read_fai_empty(self, mock_file):
        """Test reading empty FAI file."""
        mock_file.return_value.__iter__.return_value = []

        result = read_fai('empty.fai')

        assert result == {}

    @patch('builtins.open', new_callable=mock_open)
    def test_read_fai_malformed_line(self, mock_file):
        """Test reading FAI with malformed line."""
        fai_content = (
            'chr1\t1000\t5\t80\t81\nmalformed_line\nchr2\t2000\t1010\t80\t81\n'
        )
        mock_file.return_value.__iter__.return_value = fai_content.splitlines()

        result = read_fai('test.fai')

        # Should skip malformed line and continue
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

        assert result == []


class TestIsMotifInClip:
    """Test motif detection in clipped sequences."""

    def test_is_motif_in_clip_match(self):
        """Test detecting motif in clipped sequence."""
        clip_seq = 'TTAGGGTTAGGGTTAGGG'
        patterns = ['TTAGGG']

        result = isMotifInClip(clip_seq, patterns)
        assert result is True

    def test_is_motif_in_clip_no_match(self):
        """Test no motif detection in clipped sequence."""
        clip_seq = 'ATCGATCGATCG'
        patterns = ['TTAGGG']

        result = isMotifInClip(clip_seq, patterns)
        assert result is False

    def test_is_motif_in_clip_multiple_patterns(self):
        """Test motif detection with multiple patterns."""
        clip_seq = 'CCCTAACCCTAA'
        patterns = ['TTAGGG', 'CCCTAA']

        result = isMotifInClip(clip_seq, patterns)
        assert result is True  # Should match CCCTAA

    def test_is_motif_in_clip_case_insensitive(self):
        """Test case insensitive motif detection."""
        clip_seq = 'ttagggttaggg'
        patterns = ['TTAGGG']

        result = isMotifInClip(clip_seq, patterns)
        assert result is True

    def test_is_motif_in_clip_empty_sequence(self):
        """Test motif detection in empty sequence."""
        clip_seq = ''
        patterns = ['TTAGGG']

        result = isMotifInClip(clip_seq, patterns)
        assert result is False

    def test_is_motif_in_clip_empty_patterns(self):
        """Test motif detection with no patterns."""
        clip_seq = 'TTAGGGTTAGGG'
        patterns = []

        result = isMotifInClip(clip_seq, patterns)
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

    @patch('builtins.open', new_callable=mock_open)
    def test_fasta_write_read_consistency(self, mock_file):
        """Test FASTA write and read consistency."""
        # Mock file content for reading
        fasta_content = '>test_seq\nATCGATCG\n'
        mock_file.return_value.__iter__.return_value = fasta_content.splitlines()

        # Test reading
        sequences = fasta2dict('test.fasta')

        # Test writing (mock the file operations)
        mock_file.reset_mock()
        writefasta('output.fasta', 'test_seq', sequences['test_seq'])

        # Verify write was called
        mock_file.assert_called_with('output.fasta', 'w')

    def test_mask_and_filter_integration(self):
        """Test mask creation and filtering integration."""
        data = ['seq1', 'seq2', 'seq3', 'seq4', 'seq5']
        kill_indices = [1, 3]  # Remove seq2 and seq4

        # Create mask
        mask = makeMask(kill_indices, len(data))

        # Filter data
        filtered_data = filterList(data, mask)

        # Should have removed the specified indices
        expected = ['seq1', 'seq3', 'seq5']
        assert filtered_data == expected
