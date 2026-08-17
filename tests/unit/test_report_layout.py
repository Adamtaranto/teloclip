"""
Unit tests for the overhang alignment column layout.

A one-column error here still produces a plausible-looking alignment, so the
placements are pinned explicitly rather than eyeballed. The cases that matter
are the ones a flat offset gets wrong: deletions, insertions, and insertions
introduced by one read shifting every other row.
"""

from teloclip.report.layout import (
    build_columns,
    place_read,
    render_reference,
    render_row,
    terminus_column,
)

# Contig positions 1..20.
REF = 'ACGTACGTACGTACGTACGT'
LEN = 20


def lay_out(reads, is_left, window=12):
    """
    Place a set of reads and render them over shared columns.

    Parameters
    ----------
    reads : list of dict
        Read specs with cigar, seq, pos, end and clip.
    is_left : bool
        Which contig end.
    window : int, optional
        Anchored reference bases to render (default: 12).

    Returns
    -------
    dict
        Rendered reference, per-read text, columns and terminus index.
    """
    placements = [
        place_read(
            cigar=r['cigar'],
            sequence=r['seq'],
            aln_start=r['pos'],
            aln_end=r['end'],
            contig_length=LEN,
            clip_len=r['clip'],
            is_left=is_left,
            window=window,
            max_overhang=100,
        )
        for r in reads
    ]

    if is_left:
        offsets, ref_seq = list(range(0, window)), REF[:window]
    else:
        offsets, ref_seq = list(range(-window + 1, 1)), REF[-window:]

    columns = build_columns(placements, offsets)
    ref_text, ref_lead = render_reference(ref_seq, offsets, columns)

    rendered = {}
    for spec, placed in zip(reads, placements):
        lead, runs = render_row(placed, columns)
        rendered[spec['name']] = {
            'lead': lead,
            'text': ''.join(text for _, text in runs),
            'kinds': ''.join(kind[0] * len(text) for kind, text in runs),
            'runs': runs,
        }

    return {
        'columns': columns,
        'reference': ref_text,
        'reference_lead': ref_lead,
        'reads': rendered,
        'terminus': terminus_column(columns, is_left),
    }


class TestCigarAwarePlacement:
    """Indels must move a read relative to the reference, not just the clip."""

    def test_clean_alignment_tracks_the_reference(self):
        """A gapless read renders one base per reference column."""
        out = lay_out(
            [
                {
                    'name': 'clean',
                    'cigar': '5S12M',
                    'seq': 'TTTTT' + REF[:12],
                    'pos': 1,
                    'end': 12,
                    'clip': 5,
                }
            ],
            is_left=True,
        )

        # 5 clip columns outside the contig plus the 12 reference columns.
        assert len(out['columns']) == 17
        assert out['reference'] == REF[:12]
        assert out['reads']['clean']['text'] == 'TTTTT' + REF[:12]
        assert out['reads']['clean']['kinds'] == 'c' * 5 + 'a' * 12

    def test_deletion_leaves_the_read_gapped(self):
        """A deletion gaps the read while the reference keeps its bases.

        Pasting the read sequence next to the reference would shift everything
        after the deletion left by two columns.
        """
        out = lay_out(
            [
                {
                    'name': 'del2',
                    'cigar': '5S4M2D6M',
                    'seq': 'GGGGG' + REF[0:4] + REF[6:12],
                    'pos': 1,
                    'end': 12,
                    'clip': 5,
                }
            ],
            is_left=True,
        )

        read = out['reads']['del2']
        # 4 aligned, 2 deleted (gap), 6 aligned.
        assert read['text'] == 'GGGGG' + REF[0:4] + '--' + REF[6:12]
        assert read['kinds'] == 'c' * 5 + 'a' * 4 + 'g' * 2 + 'a' * 6
        # The reference is unaffected by another row's deletion.
        assert out['reference'] == REF[:12]

    def test_insertion_opens_a_column_in_every_row(self):
        """An insertion in one read gaps the reference and all other reads."""
        out = lay_out(
            [
                {
                    'name': 'plain',
                    'cigar': '5S12M',
                    'seq': 'TTTTT' + REF[:12],
                    'pos': 1,
                    'end': 12,
                    'clip': 5,
                },
                {
                    'name': 'ins3',
                    'cigar': '5S4M3I8M',
                    'seq': 'CCCCC' + REF[0:4] + 'NNN' + REF[4:12],
                    'pos': 1,
                    'end': 12,
                    'clip': 5,
                },
            ],
            is_left=True,
        )

        # 5 clip + 12 reference + 3 insertion columns.
        assert len(out['columns']) == 20
        # The inserting read fills them...
        assert out['reads']['ins3']['text'] == 'CCCCC' + REF[0:4] + 'NNN' + REF[4:12]
        # ...and everyone else, including the reference, is gapped there.
        assert out['reference'] == REF[0:4] + '---' + REF[4:12]
        assert out['reads']['plain']['text'] == 'TTTTT' + REF[0:4] + '---' + REF[4:12]

    def test_rows_stay_aligned_across_mixed_indels(self):
        """Every rendered row spans the same number of columns."""
        out = lay_out(
            [
                {
                    'name': 'clean',
                    'cigar': '5S12M',
                    'seq': 'TTTTT' + REF[:12],
                    'pos': 1,
                    'end': 12,
                    'clip': 5,
                },
                {
                    'name': 'del2',
                    'cigar': '5S4M2D6M',
                    'seq': 'GGGGG' + REF[0:4] + REF[6:12],
                    'pos': 1,
                    'end': 12,
                    'clip': 5,
                },
                {
                    'name': 'ins3',
                    'cigar': '5S4M3I8M',
                    'seq': 'CCCCC' + REF[0:4] + 'NNN' + REF[4:12],
                    'pos': 1,
                    'end': 12,
                    'clip': 5,
                },
            ],
            is_left=True,
        )

        widths = {name: r['lead'] + len(r['text']) for name, r in out['reads'].items()}
        assert len(set(widths.values())) == 1, widths
        assert out['reference_lead'] + len(out['reference']) == max(widths.values())


class TestTerminusColumn:
    """The boundary marker has to sit on the contig edge, not near it."""

    def test_left_terminus_precedes_the_first_contig_base(self):
        """At the 5' end the boundary is the left edge of offset 0."""
        out = lay_out(
            [
                {
                    'name': 'r',
                    'cigar': '5S12M',
                    'seq': 'TTTTT' + REF[:12],
                    'pos': 1,
                    'end': 12,
                    'clip': 5,
                }
            ],
            is_left=True,
        )
        # 5 clip columns come first, so the contig starts at column 5.
        assert out['terminus'] == 5

    def test_right_terminus_follows_the_last_contig_base(self):
        """At the 3' end the boundary is the right edge of offset 0."""
        out = lay_out(
            [
                {
                    'name': 'r',
                    'cigar': '12M5S',
                    'seq': REF[8:20] + 'TTTTT',
                    'pos': 9,
                    'end': 20,
                    'clip': 5,
                }
            ],
            is_left=False,
        )
        # 12 reference columns, then the clip.
        assert out['terminus'] == 12

    def test_insertion_at_the_terminus_does_not_split_the_contig(self):
        """Bases inserted at the terminal base stay on the contig side."""
        out = lay_out(
            [
                {
                    'name': 'ins',
                    'cigar': '12M2I5S',
                    'seq': REF[8:20] + 'NN' + 'TTTTT',
                    'pos': 9,
                    'end': 20,
                    'clip': 5,
                }
            ],
            is_left=False,
        )
        # 12 reference columns plus the 2 insertion columns attached to the
        # terminal base; the clip begins after all of them.
        assert out['terminus'] == 14


class TestEndsAreMirrorImages:
    """The two termini must be laid out symmetrically."""

    def test_equivalent_reads_sit_the_same_distance_outside(self):
        """A read and its mirror reach equally far past the contig."""
        left = lay_out(
            [
                {
                    'name': 'r',
                    'cigar': '6S12M',
                    'seq': 'TTTTTT' + REF[:12],
                    'pos': 1,
                    'end': 12,
                    'clip': 6,
                }
            ],
            is_left=True,
        )
        right = lay_out(
            [
                {
                    'name': 'r',
                    'cigar': '12M6S',
                    'seq': REF[8:20] + 'TTTTTT',
                    'pos': 9,
                    'end': 20,
                    'clip': 6,
                }
            ],
            is_left=False,
        )

        # Left: 6 clip columns before the terminus.
        assert left['terminus'] == 6
        assert left['reads']['r']['kinds'].startswith('c' * 6)
        # Right: 6 clip columns after it.
        assert right['terminus'] == 12
        assert right['reads']['r']['kinds'].endswith('c' * 6)


class TestSlicedReads:
    """Reads are stored truncated; offsets must still resolve."""

    def test_right_end_slice_offset_is_honoured(self):
        """A read kept only from its tail places identically to the full read."""
        full_seq = REF[8:20] + 'TTTTT'
        full = place_read(
            cigar='12M5S',
            sequence=full_seq,
            aln_start=9,
            aln_end=20,
            contig_length=LEN,
            clip_len=5,
            is_left=False,
            window=12,
            max_overhang=100,
        )
        # Same read with the first 4 bases dropped, as the capture would do.
        sliced = place_read(
            cigar='12M5S',
            sequence=full_seq[4:],
            aln_start=9,
            aln_end=20,
            contig_length=LEN,
            clip_len=5,
            is_left=False,
            window=12,
            max_overhang=100,
            seq_offset=4,
        )

        # The retained portion agrees; only the dropped prefix is missing.
        for col, base in sliced.anchor.items():
            assert full.anchor[col] == base
        assert sliced.clip == full.clip


class TestClipTruncation:
    """Long clips are trimmed toward the terminus, not away from it."""

    def test_left_clip_keeps_the_terminus_adjacent_end(self):
        """The bases nearest the contig survive truncation at the 5' end."""
        placed = place_read(
            cigar='40S12M',
            sequence='X' * 30 + 'ABCDEFGHIJ' + REF[:12],
            aln_start=1,
            aln_end=12,
            contig_length=LEN,
            clip_len=40,
            is_left=True,
            window=12,
            max_overhang=10,
        )
        clip = ''.join(placed.clip[c] for c in sorted(placed.clip))
        assert clip == 'ABCDEFGHIJ'

    def test_right_clip_keeps_the_terminus_adjacent_end(self):
        """The bases nearest the contig survive truncation at the 3' end."""
        placed = place_read(
            cigar='12M40S',
            sequence=REF[8:20] + 'ABCDEFGHIJ' + 'X' * 30,
            aln_start=9,
            aln_end=20,
            contig_length=LEN,
            clip_len=40,
            is_left=False,
            window=12,
            max_overhang=10,
        )
        clip = ''.join(placed.clip[c] for c in sorted(placed.clip))
        assert clip == 'ABCDEFGHIJ'
