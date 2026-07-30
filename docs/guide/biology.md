# Biological cases that cause trouble

Teloclip's arithmetic can be correct and its output still be wrong, because the
biology at a contig end may not be what the method assumes. This page lists the
cases that bite in practice, how to recognise them, and what to do.

The underlying assumption is: **a soft clip at a contig end represents real
sequence the assembler declined to place.** Every case below is a way that
assumption fails.

## Circular genomes

Mitochondria, chloroplasts, plasmids and nitroplasts are circular. A circular
molecule has no ends, so reads spanning the origin will *always* appear
soft-clipped at both ends of the linearised contig, at high depth, with clips
that look like perfect continuations.

They are perfect continuations — of the other end of the same molecule.
Extending them duplicates sequence.

**Recognise it:** very high overhang depth at both ends of a short contig; the
clipped sequence matches the opposite end of the same contig.

**Do:** exclude them by name.

```bash
teloclip extend reads.bam ref.fa --exclude-contigs chrM,chrC
```

Unknown names are now a hard error, so a typo will not silently do nothing.

## Collapsed repeats and rDNA arrays

Tandem arrays — rDNA especially — are usually collapsed by the assembler into
one or a few copies. Reads from every copy in the array then pile onto that
single locus. At the contig end, this looks like overwhelming overhang support.

The clipped sequence is real, but it is the *next repeat unit*, not a telomere.
Extending adds one arbitrary copy of a repeat that occurs hundreds of times.

**Recognise it:** overhang depth far above the rest of the assembly. Teloclip
flags these under **Contigs With Anomalous Overhang Coverage** in the stats
report, and logs a histogram of per-contig-end depth where they show as a long
right tail.

**Do:** inspect with `extract`, and exclude if the clipped sequence is a repeat
unit rather than telomeric.

## Interstitial telomeric sequence

Telomeric repeats occur inside chromosomes as well as at their ends —
relics of ancient fusions, or genuine internal arrays. A contig that breaks near
an interstitial telomeric array will present clipped reads full of `TTAGGG`,
which passes every motif filter.

**Recognise it:** the contig end has telomeric motifs but the contig is not a
chromosome end — check whether the other side of the break is also assembled, or
whether a linkage map or Hi-C places this end internally.

**Do:** motif presence is necessary but not sufficient. Confirm with independent
evidence before trusting an extension at a suspected internal break.

## Subtelomeric repeats

The regions just inside telomeres are often segmentally duplicated and shared
between chromosomes. A read from chromosome 3's subtelomere can align
convincingly to chromosome 7's contig end.

**Recognise it:** overhang reads at one contig end that also align elsewhere;
conflicting clipped sequences among reads at the same end.

**Do:** raise `--min-anchor` so reads need a long unique anchor, and check the
`--overhang-log` TSV for reads whose anchors are short relative to their clips.

## Homopolymer runs and low-complexity clips

ONT and PacBio both mis-call homopolymer length. A clip that is mostly a
homopolymer run carries little information and its length is unreliable.

**Recognise it:** `extend` warns when it rejects an overhang for this reason.

**Do:** `--max-homopolymer` (default 500) rejects clips containing a run at or
above that length. Lower it for noisier data.

## Chimeric reads

Long-read chimeras — two unrelated molecules ligated during library prep —
produce an alignment that stops abruptly with a large clip of unrelated
sequence. This looks exactly like a telomeric overhang except that the clipped
sequence is from somewhere else entirely.

**Recognise it:** a single read with a very long clip that disagrees with every
other read at that end.

**Do:** require motifs (`--motifs TTAGGG`), which chimeric junctions will rarely
satisfy, and prefer ends supported by several reads. Check `--min-overhangs`.

## Adapter and barcode read-through

Un-trimmed adapters appear as a consistent clipped sequence at contig ends,
identical across reads, which superficially resembles strong support.

**Recognise it:** the clipped sequence is the same short motif in every read and
matches your library's adapter.

**Do:** trim adapters before mapping. Inspect with `extract` if unsure.

## Assembly ends that are not chromosome ends

A contig end is only a candidate telomere if that contig genuinely terminates a
chromosome. Ends produced by coverage gaps, repeat collapse or misassembly will
still accumulate clipped reads.

**Do:** check the extension is consistent with your expected karyotype. A
chromosome-level assembly should have two telomeres per chromosome, not dozens
scattered across contigs. The **Contig ends extended** figure in the summary is
a quick sanity check: if it greatly exceeds twice your chromosome count, most of
those ends are not telomeres.

## A checklist before trusting an extension

- [ ] Circular contigs excluded by name
- [ ] Flagged anomalous-coverage contigs reviewed
- [ ] Clipped sequence contains the expected telomeric motif in the expected
      orientation (`CCCTAA` at the 5' end, `TTAGGG` at the 3' end for the
      canonical vertebrate repeat)
- [ ] Motif gain concentrated at contig ends, not spread through the contig
- [ ] Number of extended ends is plausible for the karyotype
- [ ] Several reads support each extended end, not one
- [ ] Extended assembly re-polished before downstream use
