## multiperipheral (in planning)

Numerically-stable 2->3 body phase space for multiperipheral (double
t-channel) diagrams, following uam19 part 5 and

    J.A.M. Vermaseren, Nucl. Phys. B229 (1983) 347-371

Motivating reactions: $e^- e^+ \to e^- e^+ \pi^0$, generalizing to
$e^- e^+ \to e^- e^+ \mu^- \mu^+$ (two-photon physics).

This directory is still in the documentation/discussion phase — no
implementation code yet. `NOTES.md` is the full running discussion log,
with every claim below anchored to quoted, line-numbered code from
`kinc2.tar.gz` (`pickin.c`, `orient.c`, `eepi.c`, `pi0.c`, `mgoto2.c`).
This file summarizes the conclusions settled so far.

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
system solve) and graph-theoretically (`networkx`). See `NOTES.md` Topic
0 for the two corrections made while getting here (a false "V>2 implies
non-tree" claim, and a `sympy` false-negative from not imposing overall
momentum conservation before checking equation redundancy).

### Internal-line momentum = signed subset-sum of external momenta

For any internal line of a tree diagram, cutting it splits the tree into
two components; conservation at every vertex on one side telescopes into
a single statement: **the momentum on that line equals (up to sign) the
sum of the external momenta on one side of the cut.** Proved rigorously
in `NOTES.md` Topic 3 from the unique-path property of trees, and
verified against our own diagram both algebraically and graph-theoretically.
This fact is the basis for the $O(2^n)$ count below — every potential
internal-line momentum in any diagram for a given process corresponds to
one subset of the external legs.

### Why potential peaks grow like $O(2^n)$ while integration variables grow like $3n-4$

For an $n$-particle final state ($n_{ext}=n+2$ external legs), the number
of *distinct possible* propagator momenta across every diagram/topology
that could contribute to the process is

$$2^{n_{ext}-1} - n_{ext} - 1 = 2^{n+1} - n - 3 \sim O(2^n)$$

(subsets of external legs, modulo the $S\leftrightarrow S^c$ symmetry,
excluding the trivial empty/full/single-leg subsets). Meanwhile the
number of non-trivial phase-space integration variables for that same
final state is $3n-4$ (Byckling–Kajantie counting: $3n$ momentum
components minus the 4 constraints of the $\delta^{(4)}$ momentum-
conservation factor; "non-trivial" because the code can additionally
avoid sampling a manifestly flat direction — e.g. `orient.c` computes the
overall-azimuth angle algebraically rather than sampling it, per Topic 2
— so `pickin`/`orient`/`eepi` consume exactly 4 random numbers for our
$n=3$ case, one under the $3n-4=5$ ceiling).

**This mismatch is a direct statement about computational cost, not just
a combinatorial curiosity** — corrected explicitly in `NOTES.md` after an
earlier draft understated this. Each of the $O(2^n)$ potential propagator
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
$O(2^n)$; it changes what has to be handled per piece.** Rather than one
universal map answerable to every potential peak, partition the
calculation (by kinematic region or by diagram) so each piece only needs
a dedicated map for its own small, bounded peak set. `pickin`/`orient`
*is* one such piece: a dedicated configuration for exactly one topology
(the double-t-channel diagram), built to tame only $t_1,t_2$ — not a
universal generator for the other 8 potential propagators that exist for
the same $n=3$ final state. The total cost of a full calculation is then
driven by the number of such pieces needed, which is itself still tied to
the diagram count — confirmed by the lecture's own aside that at 12
diagrams ($e^-e^+\to e^-e^+\mu^-\mu^+$, low energy, QED only), the
amplitude-squaring method (cost linear in diagram count) is already
preferred over generating a dedicated phase-space piece per diagram.

Full derivation, the brute-force enumeration check, and the code
evidence for "dedicated configuration" (`mgoto2.c`'s own doc comment
flagging its *lack* of protection against peaking near $\cos\theta=\pm1$,
contrasted with the lecture's description of `gamgam`'s rotate-align-
rotate-back fix) are in `NOTES.md` Topics 1 and the "Why splitting..."
section. Note: `pickin.c`'s `option` parameter is a *different*
mechanism (tunes integration order/mapping for one fixed diagram) and
must not be conflated with diagram-splitting — flagged explicitly in
`NOTES.md` to avoid repeating an earlier category of mistake.

### Figures

All diagrams in `NOTES.md` are drawn with `tikz-feynman` (LaTeX/TikZ),
following the standard particle-physics convention — straight fermion
lines, wavy photon propagators — rather than generic graph-layout tools.
Sources (`figures/*.tex`) and rendered SVGs are both version-controlled.
