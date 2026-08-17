"""
Detection of alignment file formats from their leading bytes.

Two of the sub-commands read uncompressed SAM text from a stream, and one reads
an indexed BAM. Handing either the wrong one used to fail badly: a BAM piped
into ``teloclip filter`` surfaced as a bare ``UnicodeDecodeError`` traceback
from inside the read loop, and a SAM passed to ``teloclip extend`` was reported
as a missing ``.bai``, advising the user to run ``samtools index`` on a file
that cannot be indexed.

The checks here are deliberately shallow. They read the first few bytes and
answer "is this plausibly the format we were asked for", not "is this file
valid" — pysam and the parsers remain responsible for that. The value is
entirely in the error message: every rejection names the samtools command that
would fix it.

Errors are raised as :class:`InputFormatError` rather than as click exceptions
so that this module stays free of any CLI dependency; the sub-commands
translate it into a usage error.
"""

import logging
from pathlib import Path
from typing import Iterator, Union

# gzip, and therefore BGZF, and therefore BAM. A BAM file is BGZF-compressed,
# so its first two bytes are the gzip magic number rather than the 'BAM\x01'
# that appears once the first block is decompressed.
GZIP_MAGIC = b'\x1f\x8b'

# CRAM is not compressed at the container level, so it declares itself.
CRAM_MAGIC = b'CRAM'

# How many bytes are enough to tell the formats apart. The longest magic number
# checked is four bytes; the SAM heuristic wants a little more to see a record.
_PEEK_BYTES = 64


class InputFormatError(ValueError):
    """
    Raised when an input stream or file is not in the expected format.

    The message is written for the user rather than the developer: it says what
    was found, what was expected, and the samtools invocation that converts one
    to the other.
    """


def is_gzip_magic(head: bytes) -> bool:
    """
    Report whether a byte prefix carries the gzip magic number.

    BAM is BGZF, which is gzip, so this is the test for "binary alignment
    data". It also matches a gzip-compressed SAM, which is likewise not
    readable as text and needs the same remedy.

    Parameters
    ----------
    head : bytes
        Leading bytes of the file or stream.

    Returns
    -------
    bool
        True if the prefix begins with the gzip magic number.
    """
    return head[:2] == GZIP_MAGIC


def is_cram_magic(head: bytes) -> bool:
    """
    Report whether a byte prefix declares a CRAM file.

    Parameters
    ----------
    head : bytes
        Leading bytes of the file or stream.

    Returns
    -------
    bool
        True if the prefix begins with ``CRAM``.
    """
    return head[:4] == CRAM_MAGIC


def looks_like_sam_text(head: bytes) -> bool:
    """
    Report whether a byte prefix looks like uncompressed SAM text.

    A SAM file either opens with a header line (``@HD``, ``@SQ``, ``@RG``,
    ``@PG`` or ``@CO``) or, if headerless, with a tab-delimited alignment
    record. Both are decodable ASCII, which is the cheap discriminator against
    any binary format.

    This is used to improve an error message, never to accept input, so a false
    negative costs nothing beyond a vaguer message.

    Parameters
    ----------
    head : bytes
        Leading bytes of the file or stream.

    Returns
    -------
    bool
        True if the prefix is plausibly the start of a SAM file.
    """
    if not head:
        return False

    try:
        text = head.decode('ascii')
    except UnicodeDecodeError:
        return False

    first_line = text.split('\n', 1)[0]

    if first_line.startswith(('@HD', '@SQ', '@RG', '@PG', '@CO')):
        return True

    # Headerless SAM: a QNAME, then a tab. Requiring the tab keeps plain prose
    # from being mistaken for an alignment record.
    return '\t' in first_line


def peek_bytes(stream, count: int = _PEEK_BYTES) -> bytes:
    """
    Read the leading bytes of a text stream without consuming them.

    Click hands the sub-commands a text-mode stream wrapping a buffered binary
    reader, for both real files and stdin. Peeking at that underlying buffer
    leaves the stream position untouched, so the caller can still read the
    whole file afterwards — which matters when the input is a pipe and there is
    no seeking back.

    Parameters
    ----------
    stream : file-like
        A text-mode stream, ordinarily one produced by ``click.File('r')``.
    count : int, optional
        Maximum number of bytes to look at.

    Returns
    -------
    bytes
        The leading bytes, which may be shorter than ``count`` for a short
        input, and empty if the stream exposes no peekable binary buffer (an
        in-memory ``StringIO``, for instance).
    """
    buffer = getattr(stream, 'buffer', None)
    if buffer is None:
        return b''

    try:
        return buffer.peek(count)[:count]
    except (AttributeError, OSError, ValueError):
        # Not every buffered object supports peek, and a closed or non-seekable
        # one may refuse. The check is an improvement to error reporting, so
        # declining to perform it must never itself be fatal.
        return b''


def describe_binary_input(head: bytes) -> str:
    """
    Name the binary format a byte prefix appears to be.

    Parameters
    ----------
    head : bytes
        Leading bytes of the file or stream.

    Returns
    -------
    str
        A short human-readable description of the detected format.
    """
    if is_cram_magic(head):
        return 'CRAM'
    if is_gzip_magic(head):
        return 'BAM or gzip-compressed data'
    return 'binary data'


def ensure_text_sam(stream, command: str = 'filter') -> None:
    """
    Reject a stream that is not uncompressed SAM text.

    Parameters
    ----------
    stream : file-like
        Text-mode stream about to be parsed as SAM.
    command : str, optional
        Name of the sub-command, so the suggested pipeline is the one the user
        was actually trying to run.

    Raises
    ------
    InputFormatError
        If the stream begins with gzip or CRAM magic. A stream whose leading
        bytes cannot be inspected is allowed through, and is caught later by
        :func:`sam_decode_error_message` if it turns out to be binary.
    """
    head = peek_bytes(stream)
    if not head:
        return

    if is_gzip_magic(head) or is_cram_magic(head):
        raise InputFormatError(
            f'Input to `teloclip {command}` appears to be '
            f'{describe_binary_input(head)}, but this sub-command reads '
            f'uncompressed SAM text. Convert it as you pipe:\n'
            f'    samtools view -h alignments.bam | teloclip {command} ...'
        )


def sam_decode_error_message(command: str = 'filter') -> str:
    """
    Build the message shown when SAM parsing hits undecodable bytes.

    A backstop for input that :func:`ensure_text_sam` could not inspect, such as
    a stream with no peekable buffer. The failure surfaces mid-parse as a
    decoding error, which says nothing useful on its own.

    Parameters
    ----------
    command : str, optional
        Name of the sub-command being run.

    Returns
    -------
    str
        An explanation naming the likely cause and the fix.
    """
    return (
        f'Input to `teloclip {command}` is not decodable as text, which usually '
        f'means a BAM or CRAM file was supplied where SAM was expected. '
        f'Convert it as you pipe:\n'
        f'    samtools view -h alignments.bam | teloclip {command} ...'
    )


def iter_sam_lines(stream, command: str = 'filter') -> Iterator[str]:
    """
    Iterate the lines of a SAM stream, reporting binary input clearly.

    A backstop behind :func:`ensure_text_sam`, which cannot inspect every kind
    of stream. Decoding happens lazily as lines are pulled, so binary input that
    slipped past the up-front check surfaces here, part way through parsing, as
    a decoding error that says nothing about SAM or BAM.

    Parameters
    ----------
    stream : file-like
        Text-mode stream of SAM records.
    command : str, optional
        Name of the sub-command being run, used in the error message.

    Yields
    ------
    str
        Each line of the stream, including its trailing newline.

    Raises
    ------
    InputFormatError
        If the stream cannot be decoded as text.
    """
    iterator = iter(stream)
    while True:
        try:
            line = next(iterator)
        except StopIteration:
            return
        except UnicodeDecodeError as error:
            raise InputFormatError(sam_decode_error_message(command)) from error
        yield line


def ensure_binary_bam(path: Union[str, Path]) -> None:
    """
    Reject a file that is not BGZF-compressed alignment data.

    Parameters
    ----------
    path : Union[str, Path]
        Path to the file that should be a BAM.

    Raises
    ------
    InputFormatError
        If the file begins with something other than the gzip magic number.
        Text that looks like SAM gets a message naming the conversion; anything
        else is reported as not being a BAM at all. A file that cannot be read
        here is left alone, so that the existing not-found and permission
        errors keep their own wording.
    """
    path = Path(path)

    try:
        with open(path, 'rb') as handle:
            head = handle.read(_PEEK_BYTES)
    except OSError as error:
        logging.debug(f'Could not inspect {path} for format detection: {error}')
        return

    if is_gzip_magic(head):
        return

    if looks_like_sam_text(head):
        raise InputFormatError(
            f'{path} appears to be uncompressed SAM, but `teloclip extend` '
            f'needs a sorted, indexed BAM so it can fetch alignments per '
            f'contig. Convert it:\n'
            f'    samtools view -b {path} | samtools sort -o alignments.bam\n'
            f'    samtools index alignments.bam'
        )

    if is_cram_magic(head):
        raise InputFormatError(
            f'{path} is a CRAM file, which `teloclip extend` does not read. '
            f'Convert it:\n'
            f'    samtools view -b -T reference.fasta {path} | '
            f'samtools sort -o alignments.bam\n'
            f'    samtools index alignments.bam'
        )

    raise InputFormatError(
        f'{path} is not a BAM file: it does not start with the BGZF magic '
        f'number that every BAM begins with. `teloclip extend` needs a sorted, '
        f'indexed BAM.'
    )
