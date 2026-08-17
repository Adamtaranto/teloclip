"""
Unit tests for input format detection.

The point of these checks is not to classify files correctly in the abstract —
it is to make two specific mistakes fail with a message the user can act on.
Feeding a BAM to ``teloclip filter`` used to raise a bare ``UnicodeDecodeError``
from inside the parser, and feeding a SAM to ``teloclip extend`` used to advise
running ``samtools index`` on a file that cannot be indexed. Both are common
mistakes, because ``samtools view`` output is a pipe away from either.
"""

import gzip
import io

import pytest

from teloclip.io.formats import (
    CRAM_MAGIC,
    GZIP_MAGIC,
    InputFormatError,
    describe_binary_input,
    ensure_binary_bam,
    ensure_text_sam,
    is_cram_magic,
    is_gzip_magic,
    iter_sam_lines,
    looks_like_sam_text,
    peek_bytes,
)

SAM_HEADER = b'@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:contig01\tLN:560\n'
SAM_RECORD = b'read1\t0\tcontig01\t1\t60\t100M50S\t*\t0\t0\tACGT\t*\n'


class TestMagicNumbers:
    """The byte-prefix predicates."""

    def test_gzip_magic_detected(self):
        """BGZF, and therefore BAM, starts with the gzip magic number."""
        assert is_gzip_magic(GZIP_MAGIC + b'\x08\x04rest')

    def test_gzip_magic_not_matched_by_text(self):
        """SAM text is not mistaken for compressed data."""
        assert not is_gzip_magic(SAM_HEADER)

    def test_cram_magic_detected(self):
        """CRAM declares itself in the clear."""
        assert is_cram_magic(CRAM_MAGIC + b'\x03\x00')

    @pytest.mark.parametrize('head', [b'', b'\x1f', b'C', b'CRA'])
    def test_truncated_prefixes_are_not_matched(self, head):
        """
        A prefix shorter than the magic number matches nothing.

        Parameters
        ----------
        head : bytes
            Truncated byte prefix.
        """
        assert not is_gzip_magic(head)
        assert not is_cram_magic(head)


class TestSamTextHeuristic:
    """Recognising uncompressed SAM, used only to sharpen error messages."""

    @pytest.mark.parametrize(
        'head',
        [
            SAM_HEADER,
            b'@SQ\tSN:contig01\tLN:560\n',
            b'@PG\tID:minimap2\n',
            b'@CO\ta comment\n',
            SAM_RECORD,
        ],
    )
    def test_sam_shapes_recognised(self, head):
        """
        Headered and headerless SAM are both recognised.

        Parameters
        ----------
        head : bytes
            Leading bytes of a plausible SAM file.
        """
        assert looks_like_sam_text(head)

    @pytest.mark.parametrize(
        'head',
        [
            b'',
            GZIP_MAGIC + b'\x08\x04\x00\x00',
            b'>contig01\nACGTACGT\n',  # FASTA
            b'a line of prose with no tabs\n',
        ],
    )
    def test_non_sam_rejected(self, head):
        """
        Binary data, FASTA and prose are not taken for SAM.

        Requiring a tab is what keeps the last case out: a QNAME alone is
        indistinguishable from any other word.

        Parameters
        ----------
        head : bytes
            Leading bytes that should not read as SAM.
        """
        assert not looks_like_sam_text(head)

    def test_describe_names_the_format(self):
        """The description used in error messages names what was found."""
        assert describe_binary_input(CRAM_MAGIC) == 'CRAM'
        assert 'BAM' in describe_binary_input(GZIP_MAGIC + b'\x08')
        assert describe_binary_input(b'\x00\x01\x02\x03') == 'binary data'


class TestPeekBytes:
    """Peeking must not consume, or the input would be silently truncated."""

    def test_peek_does_not_consume(self, tmp_path):
        """
        The stream still yields its full contents after a peek.

        Parameters
        ----------
        tmp_path : pathlib.Path
            Pytest temporary directory fixture.
        """
        path = tmp_path / 'alignments.sam'
        path.write_bytes(SAM_HEADER + SAM_RECORD)

        with open(path, 'r') as handle:
            head = peek_bytes(handle)
            assert head.startswith(b'@HD')
            assert handle.read().encode() == SAM_HEADER + SAM_RECORD

    def test_peek_without_binary_buffer_returns_empty(self):
        """An in-memory text stream has nothing to peek at, and that is fine."""
        assert peek_bytes(io.StringIO('@HD\tVN:1.6\n')) == b''

    def test_peek_on_empty_file(self, tmp_path):
        """
        An empty input yields an empty prefix rather than an error.

        Parameters
        ----------
        tmp_path : pathlib.Path
            Pytest temporary directory fixture.
        """
        path = tmp_path / 'empty.sam'
        path.write_bytes(b'')

        with open(path, 'r') as handle:
            assert peek_bytes(handle) == b''


class TestEnsureTextSam:
    """The guard on the SAM-reading sub-commands."""

    def test_sam_passes(self, tmp_path):
        """
        Genuine SAM text is accepted without complaint.

        Parameters
        ----------
        tmp_path : pathlib.Path
            Pytest temporary directory fixture.
        """
        path = tmp_path / 'alignments.sam'
        path.write_bytes(SAM_HEADER + SAM_RECORD)

        with open(path, 'r') as handle:
            ensure_text_sam(handle)  # must not raise

    def test_gzip_rejected_with_samtools_advice(self, tmp_path):
        """
        Compressed input is rejected, and the message names the fix.

        Parameters
        ----------
        tmp_path : pathlib.Path
            Pytest temporary directory fixture.
        """
        path = tmp_path / 'alignments.bam'
        path.write_bytes(gzip.compress(SAM_HEADER + SAM_RECORD))

        with open(path, 'r') as handle:
            with pytest.raises(InputFormatError) as excinfo:
                ensure_text_sam(handle, command='filter')

        message = str(excinfo.value)
        assert 'samtools view -h' in message
        assert 'teloclip filter' in message

    def test_command_name_appears_in_message(self, tmp_path):
        """
        The suggested pipeline names the sub-command actually being run.

        Parameters
        ----------
        tmp_path : pathlib.Path
            Pytest temporary directory fixture.
        """
        path = tmp_path / 'alignments.bam'
        path.write_bytes(gzip.compress(SAM_HEADER))

        with open(path, 'r') as handle:
            with pytest.raises(InputFormatError, match='teloclip extract'):
                ensure_text_sam(handle, command='extract')

    def test_unpeekable_stream_is_allowed_through(self):
        """
        A stream that cannot be inspected is not rejected on suspicion.

        The backstop in iter_sam_lines catches it later if it really is binary.
        """
        ensure_text_sam(io.StringIO('@HD\tVN:1.6\n'))  # must not raise


class TestIterSamLines:
    """The backstop for streams the up-front check could not inspect."""

    def test_yields_lines_unchanged(self):
        """Line content and trailing newlines survive the wrapper."""
        stream = io.StringIO('@HD\tVN:1.6\nread1\t0\tcontig01\n')

        assert list(iter_sam_lines(stream)) == [
            '@HD\tVN:1.6\n',
            'read1\t0\tcontig01\n',
        ]

    def test_empty_stream_yields_nothing(self):
        """An empty input is not an error."""
        assert list(iter_sam_lines(io.StringIO(''))) == []

    def test_decode_error_becomes_input_format_error(self, tmp_path):
        """
        Undecodable bytes are reported as a probable BAM, not as a codec fault.

        Parameters
        ----------
        tmp_path : pathlib.Path
            Pytest temporary directory fixture.
        """
        path = tmp_path / 'alignments.bam'
        path.write_bytes(gzip.compress(SAM_HEADER + SAM_RECORD))

        with open(path, 'r') as handle:
            with pytest.raises(InputFormatError) as excinfo:
                list(iter_sam_lines(handle, command='filter'))

        assert 'samtools view -h' in str(excinfo.value)


class TestEnsureBinaryBam:
    """The guard on `teloclip extend`, which needs an indexed BAM."""

    def test_bgzf_passes(self, tmp_path):
        """
        A BGZF-compressed file is accepted.

        Parameters
        ----------
        tmp_path : pathlib.Path
            Pytest temporary directory fixture.
        """
        path = tmp_path / 'alignments.bam'
        path.write_bytes(gzip.compress(SAM_HEADER))

        ensure_binary_bam(path)  # must not raise

    def test_sam_rejected_with_conversion_recipe(self, tmp_path):
        """
        A SAM is rejected with the convert-sort-index recipe, not 'index it'.

        This is the whole point of the check: `samtools index` on a SAM file
        does not work, so the previous missing-.bai message sent users down a
        dead end.

        Parameters
        ----------
        tmp_path : pathlib.Path
            Pytest temporary directory fixture.
        """
        path = tmp_path / 'alignments.sam'
        path.write_bytes(SAM_HEADER + SAM_RECORD)

        with pytest.raises(InputFormatError) as excinfo:
            ensure_binary_bam(path)

        message = str(excinfo.value)
        assert 'samtools view -b' in message
        assert 'samtools sort' in message
        assert 'samtools index alignments.bam' in message

    def test_cram_rejected_with_reference_flag(self, tmp_path):
        """
        CRAM is rejected, and the recipe includes the reference it needs.

        Parameters
        ----------
        tmp_path : pathlib.Path
            Pytest temporary directory fixture.
        """
        path = tmp_path / 'alignments.cram'
        path.write_bytes(CRAM_MAGIC + b'\x03\x00' + b'\x00' * 32)

        with pytest.raises(InputFormatError) as excinfo:
            ensure_binary_bam(path)

        assert '-T reference.fasta' in str(excinfo.value)

    def test_unrecognised_binary_rejected(self, tmp_path):
        """
        Anything else is rejected as simply not being a BAM.

        Parameters
        ----------
        tmp_path : pathlib.Path
            Pytest temporary directory fixture.
        """
        path = tmp_path / 'mystery.bin'
        path.write_bytes(b'\x00\x01\x02\x03' * 8)

        with pytest.raises(InputFormatError, match='not a BAM file'):
            ensure_binary_bam(path)

    def test_unreadable_path_is_left_to_other_checks(self, tmp_path):
        """
        A missing file is not this check's business.

        The existing not-found error is clearer than anything this could say,
        so an unreadable path passes through silently.

        Parameters
        ----------
        tmp_path : pathlib.Path
            Pytest temporary directory fixture.
        """
        ensure_binary_bam(tmp_path / 'does-not-exist.bam')  # must not raise
