# Topic 6: the muon-pair "super formulae" — $\gamma^*\gamma^*\to\mu^+\mu^-$

### Source

J.A.M. Vermaseren, *Two Photon Processes at Very High Energies*, NIKHEF-H/82-15
(July 1982) — a different, earlier paper than the 1983 Nucl. Phys. B229 one
used in Topics 0–5, but by the same author, on the same reaction family, and
using the same Levi-Civita method. Section III (pp.8–13) derives the matrix
element; Appendix C (pp.33–35) supplies the Schouten-identity reduction
formulae the derivation leans on. Local copy:
`__src__/_course_/-0-uam19/twophoton.pdf`.

### Why this is a separate, harder problem from Topic 4's $\pi^0$ case

[Topic 4](topic4-pseudoscalar-coupling.md) established that $\gamma\gamma\to\pi^0$
is the *easy* case: parity forces the coupling into one Levi-Civita term, which
is individually gauge invariant, so nothing has to cancel. $\gamma^*\gamma^*\to\mu^+\mu^-$
is the case the whole numerical-stability method actually exists for: **two**
diagrams (Fig. 3 of the paper — direct and crossed muon exchange, $k_1$- and
$k_2$-channel), each gauge invariant only *jointly* with the other, not alone.

### Step 1 — writing each individual diagram in manifestly gauge-invariant form (III.3)

Naively, each diagram's amplitude (Dirac trace between the outgoing muon
spinors, one internal muon propagator $1/(k_1^2-m^2)$ or $1/(k_2^2-m^2)$) is
*not* separately gauge invariant — replacing $e_1\to q_1$ does not kill it by
itself, because the two diagrams' non-invariant pieces cancel only against
each other. The paper's Step 1 rewrites each diagram algebraically (using the
on-shell Dirac equations $\bar u(p_6)(\slashed p_6-m)=0$ and
$(\slashed p_7+m)v(p_7)=0$, plus gamma-matrix algebra and the Schouten identity)
into a form where a genuinely single-term, manifestly-gauge-invariant piece
$\varepsilon^{q_1e_\alpha\mu\rho}\varepsilon_{q_1e_\alpha\mu\rho}$-type structure
is isolated — formula (III.3):

$$\bar u(p_6)\left\{\slashed e_1\frac{\slashed k_1+m}{k_1^2-m^2}\slashed e_2
+\slashed e_2\frac{\slashed k_2+m}{k_2^2-m^2}\slashed e_1\right\}v(p_7)
= -\varepsilon^{q_1e_\alpha\mu\rho}\varepsilon_{q_1e_\beta\mu\rho}\,
\frac{\bar u(p_6)\gamma^\beta\slashed p_6\gamma_\alpha
+\chi\gamma_\alpha\slashed p_7\gamma^\beta}{(k_1^2-m^2)(k_2^2-m^2)}\,v(p_7)$$

(the paper's exact notation for the numerator's spinor structure is denser
than reproduced here — the load-bearing point is *structural*: the sum of the
two diagrams collapses onto a single common denominator
$(k_1^2-m^2)(k_2^2-m^2)$, with the gauge-cancellation now guaranteed
automatically by construction rather than needing to be verified term by term
after the fact — "even though the gauge cancellations are hidden inside the
Levi-Civita tensors, whatever cancellation is left is due to spin structure
and therefore not serious," as the paper states directly under (III.5)).

### Step 2 — squaring and reducing via Appendix C

Squaring (III.3)/(III.4) and taking the electron/positron and muon Dirac
traces produces expressions of the shape
$\varepsilon^{q_1e_1Q\mu}\varepsilon_{q_1e_1Q\rho}\times(\text{trace stuff})$,
i.e. products of two Levi-Civita tensors sharing three contracted indices
$(q_1,e_1,Q)$ with $Q\equiv q_1+q_2$, but only two free indices left over. This
is not yet in "Gram determinant + minors" form because of the shared-index
structure — reducing it requires the **Schouten identity**

$$q_a\,\varepsilon_{bcde} = q_b\,\varepsilon_{acde}+q_c\,\varepsilon_{bade}
+q_d\,\varepsilon_{bcae}+q_e\,\varepsilon_{bcda}\tag{C1}$$

(the statement that five 4-vectors in 4D are always linearly dependent,
antisymmetrized) applied repeatedly, worked through explicitly in Appendix C
for one illustrative term (eqs. C2–C4) and used, with more steps, to derive
the two identities (C13)/(C14) needed specifically for the cross-terms $M_{12}$,
$M_{21}$ in (III.11). The paper is candid that these derivations are "far
from trivial" and were cross-checked with the symbolic-algebra program
Schoonschip [19] rather than trusted by hand alone — the same posture Topic
4/5's own verification work here has taken toward Vermaseren's other claims.

### Step 3 — the result: formula (III.5)/(III.10)/(III.11)

The final squared matrix element for two spin-$\tfrac12$ point fermions is
cast into the "$M_{ij}$" form (formula III.10):

$$|\mathcal M|^2 = \frac{W_1^{(1)}W_1^{(2)}M_{11}+W_1^{(1)}W_2^{(2)}M_{12}
+W_2^{(1)}W_1^{(2)}M_{21}+W_2^{(1)}W_2^{(2)}M_{22}}
{\left(q_1^2q_2^2\left[(q_1{\cdot}q_2)^2-(Q{\cdot}q_1)^2\right]\right)^2}$$

where $W_1,W_2$ are the (structure-function-generalizable) photon-vertex
weights, normalized so both equal $1$ for point particles, and
$M_{11},M_{12},M_{21},M_{22}$ (formula III.11) are each written purely in
terms of dot products and Levi-Civita contractions of $q_1,q_2,Q,p_1$ — no
individual term singled out as "the big one," structurally the same recipe
as `pi0.c`'s `part1..part4`, but now with four *cross* pieces ($M_{11}$ being
the "pure QED" piece with no Levi-Civita tensor left at all, $M_{12}/M_{21}$
each carrying one, $M_{22}$ carrying the most Levi-Civita structure) rather
than four minors of a single Gram determinant. The denominator's
$(q_1{\cdot}q_2)^2-(Q{\cdot}q_1)^2$ factor **looks** structurally like
`pi0.c`'s `la` term (both are quadratic brackets built from
$t_1,t_2,q_1{\cdot}q_2$), but checking this directly
(`scripts/muon_pair_denominator_check.py`) shows it is **not** the same
bracket: with $Q=q_1+q_2$, $(q_1{\cdot}q_2)^2-(Q{\cdot}q_1)^2=-t_1(t_1+2q_1{\cdot}q_2)$
exactly (a pure vector identity, no on-shell relation needed) — this is
manifestly asymmetric between $t_1$ and $t_2$ (it privileges $q_1$ through
$Q{\cdot}q_1$), whereas `pi0.c`'s `la`-bracket $t_1t_2-(q_1{\cdot}q_2)^2$ is
symmetric in $t_1,t_2$ by construction. On-shell
($q_1{\cdot}q_2=(\mu^2-t_1-t_2)/2$) the two brackets evaluate to genuinely
different polynomials in $t_1,t_2,\mu^2$ (script output: $t_1t_2-\mu^2t_1$
vs. $-\tfrac14(\mu^2-t_1-t_2)^2+t_1t_2$) — **this was an incorrect claim in
an earlier draft of this note, caught by writing and running the check
rather than asserting the match from the paper's formula shape alone.** The
two denominators are only *analogous* (both are Gram-type quadratic
invariants from the same double-t-channel kinematics), not the same
formula — the muon-pair reaction has an extra momentum ($Q=q_1+q_2$ enters
directly, unlike the $\pi^0$ case where only $t_1,t_2,q_1{\cdot}q_2$
appear), reflecting that the outgoing state at vertex $C$ is two particles,
not one.

### What the vendored code does and doesn't cover

This directory's `kinc2/` subset — `pickin.c`, `orient.c`, `pi0.c`,
`mapt.c`/`mapw.c`/`mapla.c` — implements only the $\pi^0$ ($2\to3$,
double-t-channel via $\ln(-t)$-mapping) case used in Topics 0–5. There is
**no vendored code for the muon-pair matrix element** (III.5)/(III.10)/(III.11)
in this repository, and it should not be assumed to exist.

Checked directly against the separate, un-vendored `kinc.tar.gz` archive (a
different, simpler bundle than `kinc2.tar.gz`, extracted read-only to
scratchpad for inspection, not committed here): its `meemm.c`/`meemmz.c` and
`feemm.c` files, despite the tempting filename match, are **not** an
implementation of this paper's two-photon formula either — `meemm.c` is a
plain $s$-channel/$t$-channel QED (and `meemmz.c` the QED+$Z$) matrix element
for $e^-e^+\to\mu^-\mu^+$ as an ordinary $2\to2$ process (function of a single
Mandelstam-like variable `t`, no $q_1,q_2$ two-photon structure at all), and
`feemm.c` drives it through `mgoto2.c` — a generic $2\to2$ two-body-decay
phase-space generator (samples $\cos\theta,\phi$ in the CM frame via
`mgoto2`'s own random-number-to-angle map), not the double-t-channel
`pickin`/`orient`/$\ln(-t_1,-t_2)$ machinery this repository's `kinc2/`
vendors. `mgoto3.c` (a three-body Dalitz-variable generator) is likewise
unrelated — used elsewhere in that archive's suite for different final
states, not for this reaction. **This was checked, not assumed**: grepped
the whole `kinc.tar.gz` archive for `part1`, `part2`, `dd4`, `dd2` (the
naming pattern `pi0.c` uses for Gram-determinant/minor pieces) and found no
matches anywhere — confirming this simpler archive genuinely does not carry
an implementation of formula (III.5)/(III.10)/(III.11).

**So this topic's derivation is grounded directly in the paper's own
formulae and its Appendix C reduction identities, not cross-checked against
a driver the way Topics 4/5 cross-checked `pi0.c`.** That is an honest
difference in verification depth from the earlier topics, not glossed over:
Topics 4/5's Gram-minor and $s$-scaling claims were checked both
symbolically *and* numerically against running vendored code; this topic's
$M_{11}$–$M_{22}$ structure has been checked only symbolically (the
denominator identity above) and by matching the paper's own internal
consistency statements (III.9's $W_1=W_2=1$ point-particle normalization,
III.11's dimensional/index structure). No numerical driver for this reaction
exists in the repository to check the $M_{ij}$ formulae against directly.

### Open threads / not yet resolved

- The $M_{11}$–$M_{22}$ formulae (III.11) have not been checked numerically
  against a running implementation — none exists in this repository or the
  vendored archives — only symbolically, and only for the denominator
  structure, not the full numerator reduction chain (C1)→(C13)/(C14).
- Appendix C's derivation chain (C2)→(C14) has been read and structurally
  summarized here, but not independently re-derived term-by-term the way
  Topic 4's Gram-minor claim was (there, a from-scratch symbolic
  verification in `scripts/gram_minors_pi0.py` reproduced the claim
  independently of the paper's own algebra; here, the reduction identities
  are taken on the paper's authority, cross-checked only for
  self-consistency of the final formula, not re-derived from the Schouten
  identity from scratch).
- Writing an actual `pickin`/`orient`-style driver for this reaction (a
  genuine $2\to4$ process, $e^-e^+\to e^-e^+\mu^-\mu^+$, one order more
  complex than the $\pi^0$'s $2\to3$) was explicitly out of scope here and
  remains undone — it is the natural "close the loop" step if this thread is
  picked up again, mirroring what `check_gram_minors.c` did for Topic 4.
