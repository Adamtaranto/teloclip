# Why coverage dies at contig ends

Teloclip exists because of a sampling artefact: read coverage does not stay
flat all the way to the end of a linear chromosome — it decays smoothly to
zero over roughly one read length. Assemblers see thin, one-sided evidence in
that window and stop early, which is why draft contigs so often end just short
of the telomere. This page gives the statistical model behind that decay, so
you can reason about how much overhang evidence to expect at a given depth
and read-length profile.

If you have not already, read [how overhangs are defined](overhangs.md) and
[telomere biology](biology.md) first — this page explains why the overhang
reads those pages describe are always rarer than the interior depth suggests.

## The model

Two random processes place a read on a molecule:

1. **Fragmentation is uniform.** Breakpoints fall anywhere along the
   molecule, so a fragment of length \(\ell\) can start at any position that
   lets it fit.
2. **Size selection is a hard cutoff.** Library preparation (or a read-length
   filter at mapping time) discards fragments shorter than some minimum
   \(L_{\min}\).

Consider a base at distance \(x\) from the terminus (\(x = 0\) is the
outermost base). A fragment of length \(\ell\) covers it only if the fragment
starts within the last \(x\) bases before it — there are
\(\min(x + 1, \ell)\) such placements, against \(\ell\) placements for an
interior base. A terminal base is simply a smaller target.

Averaging over the fragment-length distribution (truncated at
\(L_{\min}\)) gives the expected coverage relative to the interior:

$$
c(x) \;=\; \frac{\mathbb{E}\left[\min(x + 1,\, \ell) \mid \ell \ge L_{\min}\right]}
                 {\mathbb{E}\left[\ell \mid \ell \ge L_{\min}\right]}
$$

Note the denominator: size selection also changes the *interior* mean
fragment length, so both numerator and denominator use the truncated
distribution. Absolute expected coverage is \(D \cdot c(x)\), where \(D\) is
the interior depth after size selection.

The model assumes the molecule is far longer than any fragment
(semi-infinite): on a real chromosome the same decay is simply mirrored at
the other terminus.

### Fixed fragment length

If every fragment has length \(L\), the curve is a straight ramp:

$$
c(x) = \frac{\min(x + 1,\, L)}{L}
$$

Coverage rises linearly from \(1/L\) at the outermost base and reaches the
interior depth exactly one fragment length in.

### Lognormal fragment lengths

Real long-read length profiles are well described by a lognormal. With
\(\ell \sim \mathrm{Lognormal}(\mu, \sigma)\) (parameterised by the
real-space mean \(m\) and standard deviation \(s\)), the partial expectation
has a closed form in the normal CDF \(\Phi\):

$$
\mathbb{E}\left[\min(a, \ell)\right]
  = m\,\Phi\!\left(\frac{\ln a - \mu - \sigma^2}{\sigma}\right)
  + a\left(1 - \Phi\!\left(\frac{\ln a - \mu}{\sigma}\right)\right)
$$

so \(c(x)\) is cheap to evaluate everywhere — no simulation needed. The decay
is no longer linear: the long tail of the distribution reaches deeper into
the contig, while short fragments (if they survive selection) soften the
curve near the terminus.

## What the curves look like

The figures below use a PacBio HiFi-like library: lognormal fragment lengths
with mean 15 kb and sd 3 kb, and a strict lower-bound size-selection cutoff
of 10 kb (every fragment shorter than the cutoff is discarded).

![Two-panel figure of analytic coverage-decay curves. Panel a: relative coverage against distance from the contig end for a lognormal read profile with mean 15 kb and sd 3 kb under lower-bound size-selection cutoffs of none, 10, 15 and 20 kb; harsher cutoffs push the decay region deeper into the contig. Panel b: the same curve at the 10 kb cutoff for the HiFi profile, a wider sd 6 kb lognormal, and a fixed 15 kb fragment length; wider length distributions stretch the decay further.](../images/coverage-decay-analytic.png)

Three things to take away:

- **The decay region is about one (truncated) mean fragment length deep.**
  With the 15 ± 3 kb HiFi profile, coverage within 5 kb of the terminus is
  only about a third of the interior depth, and interior depth is not
  reached until ~20 kb in.
- **Below the cutoff distance the curve is an exact straight line.** Every
  surviving fragment is at least 10 kb, so for \(x < L_{\min}\) the number
  of placements that cover \(x\) grows linearly with \(x\) regardless of the
  length distribution — the distribution only shapes the shoulder beyond the
  cutoff.
- **Cutoffs bite when they reach the typical length.** The 10 kb HiFi cutoff
  discards only ~3% of a 15 ± 3 kb library, so its curve sits almost on the
  no-cutoff curve; pushing the cutoff to 15–20 kb (at or above the mean)
  discards most fragments and drags the depleted region far deeper.

Simulation confirms the analytic curve. Each simulation fragments a 1 Mb
molecule with lognormal fragment lengths (mean 15 kb, sd 3 kb), discards
fragments below the 10 kb lower-bound cutoff, places the survivors at 30×
interior depth, and tallies terminal coverage; the band spans the central
95% of 200 independent runs:

![Simulated terminal coverage for lognormal fragment lengths with mean 15 kb and standard deviation 3 kb under a strict 10 kb lower-bound size-selection cutoff at thirty-fold interior depth, with a 95 percent bootstrap band from 200 simulations, overlaid with the analytic expectation, which tracks the simulated mean closely across the whole decay region.](../images/coverage-decay-bootstrap.png)

!!! note "Reproducing these figures"

    Both figures are generated by
    [`scripts/generate_model_figures.py`](https://github.com/adamtaranto/teloclip/blob/main/scripts/generate_model_figures.py)
    using `teloclip.model` (install with `pip install 'teloclip[model]'`).
    The analytic curve, the simulator and the bootstrap band are all
    available as library functions in
    [`teloclip.model.decay`](https://github.com/adamtaranto/teloclip/blob/main/src/teloclip/model/decay.py).

## Explore the parameters

Use the sliders to see how depth, read-length profile and size-selection
cutoff shape the decay. The curve is the analytic \(c(x)\) computed in your
browser; the readout reports expected absolute coverage near the terminus.

<div id="coverage-decay-app" markdown="0">
  <noscript>Interactive figure requires JavaScript; the static figures above
  show the same model.</noscript>
</div>

## What this means for teloclip

- **Overhang reads are the survivors of this decay.** A read can only carry
  an overhang if it covers the outermost base, where coverage is at its
  minimum — at 30× interior depth with the 15 ± 3 kb HiFi profile and a
  10 kb cutoff, expect only a handful of reads to span any given terminus. Zero overhangs
  at one contig end is often just sampling, not a structural problem.
- **Depth for telomere recovery is set by \(D \cdot c(0^+)\), not \(D\).**
  If you want, say, ≥5 candidate overhang reads per end, work back from the
  curve: interior depth must be several times higher than a naive estimate.
- **Be gentle with length filters.** Mapping-time read-length cutoffs act
  exactly like library size selection in this model: every kilobase you add
  to the minimum read length pushes the depleted zone deeper and removes
  potential overhang evidence. Filter on alignment quality and anchor length
  (see [how overhangs are defined](overhangs.md)) rather than raw length
  where you can.
- **Interpret the report's depth panels with the curve in mind.** The
  [extend report](reports.md) plots overhang depth per contig end; ends with
  few overhangs but healthy flanking coverage are consistent with this model,
  while ends where coverage collapses *faster* than the curve predicts may
  indicate misassembly or collapsed repeats.
