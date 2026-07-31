# How overhangs are defined

Every teloclip sub-command asks the same question of every alignment: *is this
soft clip a valid overhang at a contig end?* This page defines that precisely,
because the answer determines which reads survive filtering and which sequence
ends up in your assembly.

All three sub-commands share one implementation, so `filter`, `extract` and
`extend` agree by construction.

!!! note "This changed"

    Before v0.3.5 each sub-command had its own copy of these rules, and the
    copies disagreed by a base at each end and about whether `--min-clip`
    applied at all. See the [CHANGELOG](https://github.com/adamtaranto/teloclip/blob/main/CHANGELOG.md)
    if you are comparing results across versions.

## Coordinates

Reference coordinates are **1-based and inclusive**, matching the SAM `POS`
field.

```
ref_span  = sum of CIGAR lengths over M, D, N, =, X
aln_start = POS                     # first aligned reference base
aln_end   = POS + ref_span - 1      # last aligned reference base
```

`aln_end` is the last aligned base, not one position past it. Insertions (`I`)
and soft clips (`S`) consume read bases but not reference bases, so they do not
contribute to the span.

## Gap and overhang

For a contig of length `L`:

```
gap_left  = aln_start - 1     # terminal contig bases the alignment does not cover
gap_right = L - aln_end

overhang  = clip_len - gap    # clipped read bases lying past the contig end
```

The soft clip is contiguous with the alignment in read space, so its bases
notionally occupy the reference positions immediately before `aln_start` (or
after `aln_end`). Some of those fall *inside* the contig — that is the `gap` —
and the rest fall outside it.

```
                gap = 3        overhang = 5
             |<------->|<--------------->|
contig       |=========|===================================>
read              ~~~~~~~~~~~~~~~----------------------->
             clip_len = 8      aln_start
```

## Acceptance

Applied identically at both ends:

| Criterion | Rule |
|---|---|
| `--max-break` | `gap <= max_break` — **inclusive** |
| `--min-clip` | `overhang >= min_clip` |
| `--min-anchor` | `anchor >= min_anchor`, where anchor counts `M`/`=`/`X` bases only |

Two consequences worth internalising:

**`--max-break 0` means flush.** The alignment must reach the terminal base
exactly: `aln_start == 1` on the left, `aln_end == L` on the right.

**`--min-clip` measures novel sequence, not clip length.** A 500 bp clip sitting
499 bp inside the contig contributes exactly 1 bp of new sequence. Thresholding
the raw clip length would treat it as a strong candidate; thresholding the
overhang treats it as the marginal one it is. This matches the option's
long-standing help text.

## Why this guarantees contigs never shrink

Applying an overhang trims the `gap` bases the read did not cover, then grafts
on the **whole** clip in their place. The trimmed region is replaced by the
read's own version of it, plus the novel tail. So:

```
net change in contig length = overhang
```

Because `--min-clip` is clamped to at least 1, every accepted overhang has
`overhang >= 1`, and an extension that would shorten a contig cannot be
represented. `extend` additionally re-checks the net gain immediately before
writing, as a backstop.

## Anchoring

The anchor is the number of read bases actually aligned to the reference:
`M`, `=` and `X` only.

Deletions advance the reference but contribute no read bases; insertions
contribute read bases with no reference support. Neither is evidence that the
read is correctly placed, so neither counts. A read whose "alignment" is mostly
insertions is not well anchored no matter how long it is.

Raise `--min-anchor` if you are seeing spurious overhangs from repetitive
regions; a read needs a substantial unique anchor to be trusted at a contig end.
