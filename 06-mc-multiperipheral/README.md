## multiperipheral (in planning)

Numerically-stable 2->3 body phase space for multiperipheral (double
t-channel) diagrams, following uam19 part 5 and

    J.A.M. Vermaseren, Nucl. Phys. B229 (1983) 347-371

Motivating reactions: $e^- e^+ \to e^- e^+ \pi^0$, generalizing to
$e^- e^+ \to e^- e^+ \mu^- \mu^+$ (two-photon physics).

This directory is still in the documentation/discussion phase for the
double-t-channel phase space itself — no standalone event generator has
been written yet. What it does contain: [`NOTES.md`](NOTES.md), an index
into the full discussion split by topic under [`notes/`](notes/), with
every claim anchored to quoted, line-numbered code from `kinc2.tar.gz`
(`pickin.c`, `orient.c`, `eepi.c`, `pi0.c`, `mgoto2.c`) and, where
relevant, to the primary sources the lecture cites (Vermaseren 1983;
Byckling & Kajantie, *Particle Kinematics*, 1973; Byers & Yang,
*Rev. Mod. Phys.* 36 (1964) 595); [`scripts/`](scripts/), standalone
Python checks of specific combinatorial/graph-theoretic claims;
[`kinc2/`](kinc2/), the subset of the vendored reference implementation
needed to build the one non-Python check so far,
[`kinc2-driver/scan_pi0_scaling.c`](kinc2-driver/scan_pi0_scaling.c), a
driver that instruments the real `pickin`/`orient`/`pi0` code directly
(see Topic 5 below). This file summarizes the conclusions settled so
far.

### Tree diagrams and loops, defined rigorously

"Tree diagram" is not "no cycle when you draw it" — that's a consequence,
not the definition. The rigorous definition is by **constraint counting**:
at each of $V$ vertices, 4-momentum conservation gives one equation
relating the lines meeting there; exactly one of those $V$ equations is
always redundant (implied by the other $V-1$, via overall momentum
conservation), so there are $V-1$ independent constraints. With $E$
internal-line momenta as unknowns,

$$L = E - (V-1)$$

is the number of internal momenta left undetermined — the number of loop
integration variables. $L=0$ is a tree; $L\geq1$ is a loop diagram. This
was checked against three worked examples (our own double-t-channel
diagram, $L=0$; a one-loop QED box, $L=1$; a 4-vertex chain diagram with
the *same* vertex count as the box but $L=0$ — the direct counterexample
to "just count vertices"), verified both symbolically (`sympy` vertex
system solve) and graph-theoretically (`networkx`). See
[`notes/topic0-tree-vs-loop.md`](notes/topic0-tree-vs-loop.md) for the two
corrections made while getting here (a false "V>2 implies non-tree"
claim, and a `sympy` false-negative from not imposing overall momentum
conservation before checking equation redundancy).

### Internal-line momentum = signed subset-sum of external momenta

For any internal line of a tree diagram, cutting it splits the tree into
two components; conservation at every vertex on one side telescopes into
a single statement: **the momentum on that line equals (up to sign) the
sum of the external momenta on one side of the cut.** Proved rigorously
in [`notes/topic3-subset-sum-proof.md`](notes/topic3-subset-sum-proof.md)
from the unique-path property of trees, and verified against our own
diagram both algebraically and graph-theoretically. This fact is the
basis for the $O(2^n)$ count below — every potential internal-line
momentum in any diagram for a given process corresponds to one subset of
the external legs. The same telescoping argument generalizes to loop
diagrams (internal-line momentum = external subset-sum + integer
combination of loop momenta), verified on a one-loop QED box example — see
[`notes/topic1-peaks-and-O2n.md`](notes/topic1-peaks-and-O2n.md)'s
"Non-tree generalization" section.

### Why potential peaks grow like $O(2^n)$ while integration variables grow like $3n-4$

For an $n$-particle final state ($n_{ext}=n+2$ external legs), the number
of *distinct possible* propagator momenta across every **tree** diagram
that could contribute to the process (this restriction to trees is
explicit and load-bearing — Vermaseren's own sentence is unrestricted, but
the proof behind this count only holds for trees; see
[`notes/topic1-peaks-and-O2n.md`](notes/topic1-peaks-and-O2n.md) for the
scope note and the loop-diagram generalization, which shows a single loop
diagram's genuine singularity count stays the same order, $O(2^E)$, but
leaves the harder question — summing over all diagram topologies at all
loop orders — open) is

$$2^{n_{ext}-1} - n_{ext} - 1 = 2^{n+1} - n - 3 \sim O(2^n)$$

(subsets of external legs, modulo the $S\leftrightarrow S^c$ symmetry,
excluding the trivial empty/full/single-leg subsets; checked to be exactly
the *achievable* bipartitions, not an overcount, by enumerating every
cubic tree topology for $n=3$). Meanwhile the number of non-trivial
phase-space integration variables for that same final state is $3n-4$ —
checked directly against Byckling & Kajantie's Table III.1 (1973, p.49),
which distinguishes four different variable counts for a $2\to n$
process; the one relevant here ("Final State Variables," $3n-4$, which
*includes* one trivial overall azimuth $\phi$) reduces to "Essential Final
State Variables" ($3n-5$) once $\phi$ is dropped — for our $n=3$ case,
$3n-5=4$, matching exactly that `pickin`/`orient`/`eepi` consume exactly 4
random numbers (`orient.c` computes the azimuth algebraically rather than
sampling it).

**This mismatch is a direct statement about computational cost, not just
a combinatorial curiosity.** Each of the $O(2^n)$ potential propagator
momenta is a place the integrand can develop an uncontrolled peak; taming
one requires a dedicated importance-sampling map (`mapt.c`'s `dt/t`,
`mapw.c`'s `dw^2/w^2`, etc.), and each such map consumes one of the
$3n-4$ available coordinate "slots." An exponentially-growing hazard
count against a linearly-growing map budget means no single
$3n-4$-dimensional coordinate system can, in general, be built tangent to
every potential peak at once — this is exactly why Vermaseren opens the
lecture with "this does not mean... we can calculate all reactions
trivially."

**The fix — splitting phase space into pieces — does not shrink
$O(2^n)$; it changes what has to be handled per piece.** Checked directly
against Vermaseren's own 1983 paper (p.349), which distinguishes two
methods: "event generation of the second kind" (rewrite variables for
experimental cuts, then split phase space into pieces such that **each
piece contains a single peak** — the paper's own words) versus "event
generation of the third kind" — **the method `pickin`/`orient` actually
implements**: use the propagator denominators $t_1,t_2$ directly as
integration variables, mapped via $\ln(-t_1),\ln(-t_2)$, handling both
peaks *in one unsplit calculation*. So `pickin`/`orient` is not an example
of splitting — it's a worked example of the paper's alternative,
universality-preferring method, used here because this one diagram only
has two peaks to tame. The `mgoto2`/`gamgam` rotate-align-rotate-back
trick described in the lecture is a third, distinct technique (a
coordinate transformation applied before sampling, not phase-space
splitting or the $\ln(-t)$ method). Full derivation, the code evidence,
and the correction of an earlier mislabeling of `pickin`/`orient` are in
[`notes/topic1-peaks-and-O2n.md`](notes/topic1-peaks-and-O2n.md). Note:
`pickin.c`'s `option` parameter is a *different* mechanism still (tunes
integration order/mapping for one fixed diagram within the $\ln(-t)$
method) and must not be conflated with any of the above.

### Why the pseudoscalar ($\pi^0$) coupling is the cleanest worked example

Choosing $e^-e^+\to e^-e^+\pi^0$ as the first worked matrix element isn't
about kinematics (identical to the QED case above) — it's the simplest
possible test case for the paper's numerical-stability *method*. The
$\pi^0\gamma\gamma$ coupling is forced by parity into a single Levi-Civita
term, which is automatically, individually gauge invariant — so every term
in the squared matrix element (`pi0.c`'s four `part1..part4`, matching the
paper's formula 3.2) is a genuine piece of the Gram-determinant family,
with nothing to cancel. The genuinely hard case — QED
$\gamma^*\gamma^*\to\mu^+\mu^-$ — needs two diagrams that are gauge
invariant only jointly, which is the actual source of the severe numerical
cancellations the whole paper exists to solve. See
[`notes/topic4-pseudoscalar-coupling.md`](notes/topic4-pseudoscalar-coupling.md).

### Why "grows with powers of $s$" is a catastrophe, and where it doesn't cancel

The lecture note warns that formula (3.2)'s first term, combined with
$1/(t_1t_2)^2$, "would give a crosssection that grows with powers of
$s$" — this is a unitarity statement: the near-real-photon corner
($t_1,t_2\to0$) is the equivalent-photon-approximation regime, where the
process factorizes into a photon flux times the fixed,
$s$-independent $\sigma(\gamma\gamma\to\pi^0)$; a term surviving with
power-law $s$-growth there would contradict that factorization. **This
does not happen via cancellation between `pi0.c`'s four terms.** An
earlier draft of Topic 4 assumed it did; checked directly (both
analytically,
[`scripts/gram_s_scaling.py`](scripts/gram_s_scaling.py), which derives
$\varepsilon^{p_1q_1p_2q_2}\varepsilon_{p_1q_1p_2q_2}\sim s^2t_1t_2$ at
fixed $t_1,t_2$, and numerically, against the real `pickin`/`orient`/`pi0`
code,
[`kinc2-driver/scan_pi0_scaling.c`](kinc2-driver/scan_pi0_scaling.c)) it
is now confirmed the sum `part1+part2+part3+part4` genuinely grows like
$s^2$ at fixed small $t_1,t_2$, tracking `part1` with no cancellation —
exactly matching the paper's own statement (p.353) that "all 4 terms are
positive so no cancellations occur." The "decent, $1/(t_1t_2)$" behavior
the lecture note promises is a property of the full cross-section
(matrix element $\times$ phase-space Jacobian $\times$ the $s$-dependent
coupling between $t_1,t_2$'s allowed range and $s$ at the phase-space
boundary), not of `pi0.c`'s bare squared-matrix-element output taken
alone. See
[`notes/topic5-s-scaling-pi0.md`](notes/topic5-s-scaling-pi0.md).

### Figures

All diagrams in the `notes/` files are drawn with `tikz-feynman`
(LaTeX/TikZ), following the standard particle-physics convention —
straight fermion lines, wavy photon propagators — rather than generic
graph-layout tools. Sources (`figures/*.tex`) and rendered SVGs are both
version-controlled.
