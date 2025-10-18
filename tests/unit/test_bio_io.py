"""
Unit tests for bio_io module.

Tests the BioPython-based memory-efficient FASTA I/O functions.
"""

import tempfile
from pathlib import Path
from unittest import mock
from unittest.mock import patch

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from teloclip.bio_io import (
    load_fasta_sequences,
    stream_fasta_sequences,
    stream_write_fasta_sequences,
    validate_fasta_against_fai,
    write_fasta_sequences,
)


class TestLoadFastaSequences:
    """Test load_fasta_sequences function."""

    def test_load_fasta_sequences_simple(self):
        """Test loading simple FASTA file."""
        fasta_content = ">seq1 description 1\nATCG\n>seq2 description 2\nGCTA\n"
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fna') as f:
            f.write(fasta_content)
            f.flush()
            
            sequences = load_fasta_sequences(f.name)
            
        assert len(sequences) == 2
        assert 'seq1' in sequences
        assert 'seq2' in sequences
        assert sequences['seq1'] == ('seq1 description 1', 'ATCG')
        assert sequences['seq2'] == ('seq2 description 2', 'GCTA')

    def test_load_fasta_sequences_multiline(self):
        """Test loading FASTA with multiline sequences."""
        fasta_content = ">seq1\nATCG\nGCTA\n>seq2\nTTTT\nAAAA\n"
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fna') as f:
            f.write(fasta_content)
            f.flush()
            
            sequences = load_fasta_sequences(f.name)
            
        assert sequences['seq1'] == ('seq1', 'ATCGGCTA')
        assert sequences['seq2'] == ('seq2', 'TTTTAAAA')

    def test_load_fasta_sequences_empty(self):
        """Test loading empty FASTA file."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fna') as f:
            f.write("")
            f.flush()
            
            sequences = load_fasta_sequences(f.name)
            
        assert len(sequences) == 0


class TestStreamFastaSequences:
    """Test stream_fasta_sequences function."""

    def test_stream_fasta_sequences_simple(self):
        """Test streaming simple FASTA file."""
        fasta_content = ">seq1 description 1\nATCG\n>seq2 description 2\nGCTA\n"
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fna') as f:
            f.write(fasta_content)
            f.flush()
            
            sequences = list(stream_fasta_sequences(f.name))
            
        assert len(sequences) == 2
        assert sequences[0] == ('seq1', 'seq1 description 1', 'ATCG')
        assert sequences[1] == ('seq2', 'seq2 description 2', 'GCTA')

    def test_stream_fasta_sequences_iterator(self):
        """Test that streaming returns an iterator, not loaded data."""
        fasta_content = ">seq1\nATCG\n>seq2\nGCTA\n"
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fna') as f:
            f.write(fasta_content)
            f.flush()
            
            seq_iter = stream_fasta_sequences(f.name)
            
            # Should be able to iterate one at a time
            first = next(seq_iter)
            assert first == ('seq1', 'seq1', 'ATCG')
            
            second = next(seq_iter)
            assert second == ('seq2', 'seq2', 'GCTA')
            
            # Should be exhausted
            with pytest.raises(StopIteration):
                next(seq_iter)


class TestWriteFastaSequences:
    """Test write_fasta_sequences function."""

    def test_write_fasta_sequences_to_file(self):
        """Test writing sequences to a file."""
        sequences = {
            'seq1': ('seq1 description 1', 'ATCG'),
            'seq2': ('seq2 description 2', 'GCTA')
        }
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fna') as f:
            output_path = f.name
            
        write_fasta_sequences(sequences, output_path)
        
        # Read back and verify
        with open(output_path, 'r') as f:
            content = f.read()
            
        assert '>seq1 description 1' in content
        assert '>seq2 description 2' in content
        assert 'ATCG' in content
        assert 'GCTA' in content

    @patch('sys.stdout')
    def test_write_fasta_sequences_to_stdout(self, mock_stdout):
        """Test writing sequences to stdout."""
        sequences = {
            'seq1': ('seq1 description', 'ATCG')
        }
        
        write_fasta_sequences(sequences, None)
        
        # Should have called SeqIO.write with stdout
        mock_stdout.write.assert_called()


class TestStreamWriteFastaSequences:
    """Test stream_write_fasta_sequences function."""

    def test_stream_write_fasta_sequences_to_file(self):
        """Test streaming write to file."""
        sequences = [
            ('seq1', 'seq1 description 1', 'ATCG'),
            ('seq2', 'seq2 description 2', 'GCTA')
        ]
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fna') as f:
            output_path = f.name
            
        stream_write_fasta_sequences(iter(sequences), output_path)
        
        # Read back and verify
        with open(output_path, 'r') as f:
            content = f.read()
            
        assert '>seq1 description 1' in content
        assert '>seq2 description 2' in content
        assert 'ATCG' in content
        assert 'GCTA' in content

    @patch('sys.stdout')
    def test_stream_write_fasta_sequences_to_stdout(self, mock_stdout):
        """Test streaming write to stdout."""
        sequences = [('seq1', 'seq1 description', 'ATCG')]
        
        stream_write_fasta_sequences(iter(sequences), None)
        
        # Should have called stdout
        mock_stdout.write.assert_called()


class TestValidateFastaAgainstFai:
    """Test validate_fasta_against_fai function."""

    def test_validate_fasta_against_fai_perfect_match(self):
        """Test validation with perfect match."""
        fasta_content = ">seq1\nATCG\n>seq2\nGCTA\n"
        fai_dict = {'seq1': 4, 'seq2': 4}
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fna') as f:
            f.write(fasta_content)
            f.flush()
            
            missing_from_fasta, missing_from_fai = validate_fasta_against_fai(f.name, fai_dict)
            
        assert len(missing_from_fasta) == 0
        assert len(missing_from_fai) == 0

    def test_validate_fasta_against_fai_missing_from_fasta(self):
        """Test validation with sequences missing from FASTA."""
        fasta_content = ">seq1\nATCG\n"
        fai_dict = {'seq1': 4, 'seq2': 4, 'seq3': 4}
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fna') as f:
            f.write(fasta_content)
            f.flush()
            
            missing_from_fasta, missing_from_fai = validate_fasta_against_fai(f.name, fai_dict)
            
        assert missing_from_fasta == {'seq2', 'seq3'}
        assert len(missing_from_fai) == 0

    def test_validate_fasta_against_fai_missing_from_fai(self):
        """Test validation with sequences missing from FAI."""
        fasta_content = ">seq1\nATCG\n>seq2\nGCTA\n>seq3\nTTTT\n"
        fai_dict = {'seq1': 4}
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fna') as f:
            f.write(fasta_content)
            f.flush()
            
            missing_from_fasta, missing_from_fai = validate_fasta_against_fai(f.name, fai_dict)
            
        assert len(missing_from_fasta) == 0
        assert missing_from_fai == {'seq2', 'seq3'}


class TestBioIoIntegration:
    """Integration tests for bio_io module."""

    def test_round_trip_consistency(self):
        """Test that load -> write -> load produces same result."""
        original_sequences = {
            'seq1': ('seq1 test sequence', 'ATCGATCGATCG'),
            'seq2': ('seq2 another sequence', 'GCTAGCTAGCTA')
        }
        
        # Write to temporary file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fna') as f:
            temp_path = f.name
            
        write_fasta_sequences(original_sequences, temp_path)
        
        # Load back
        loaded_sequences = load_fasta_sequences(temp_path)
        
        # Should be identical
        assert loaded_sequences == original_sequences

    def test_streaming_vs_loading_consistency(self):
        """Test that streaming and loading produce equivalent results."""
        fasta_content = ">seq1 description\nATCGATCG\n>seq2 other\nGCTAGCTA\n"
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fna') as f:
            f.write(fasta_content)
            f.flush()
            
            # Load all at once
            loaded = load_fasta_sequences(f.name)
            
            # Stream and collect
            streamed = {seq_id: (desc, seq) for seq_id, desc, seq in stream_fasta_sequences(f.name)}
            
        assert loaded == streamed
