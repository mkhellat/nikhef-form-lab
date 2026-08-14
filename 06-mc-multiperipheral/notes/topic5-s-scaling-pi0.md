# Topic 5: what "would grow with powers of $s$" actually refers to in formula (3.2)

### The source passage

> "The most important part is that the first term combines with
> $1/(t_1t_2)^2$ which by itself would cause problems with physics because
> it would give a crosssection that grows with powers of $s$. Only a
> behaviour like $1/(t_1t_2)$ could give a decent crosssection. This means
> that when $t_1$ and $t_2$ become small, the numerator must become very
> small as well. But the numerator, if we would contract the Levi-Civita
> tensors and write them in terms of dotproducts, has terms that are very
> big ($O(s^2M_\pi^4)$). This indicates that there must be very bad
> cancellations between the various terms."
> — uam19 `part5.pdf`, discussing Vermaseren 1983 formula (3.2)

### The question

Two things needed pinning down precisely, because an imprecise reading of
this passage leads to a wrong physical picture:

1. **Why is "a cross-section that grows with powers of $s$" a catastrophe
   as a matter of principle** — not just an inconvenient number?
2. **Where, mechanically, does the cancellation happen** that turns a
   naively $s^2$-sized numerator into the "decent," $s$-independent
   $1/(t_1t_2)$ behavior the passage says is required? An earlier draft of
   this discussion assumed the four terms of formula (3.2)
   (`pi0.c`'s `part1..part4`) cancel against each other. **That assumption
   was checked here and found wrong** — see below.

### Part 1: why $s$-growth is a catastrophe in principle

This is a statement about unitarity, not aesthetics. The full QED S-matrix
is unitary ($S^\dagger S = 1$); this bounds every partial-wave amplitude
for any fixed physical process, and hence bounds how the corresponding
cross-section can grow with $s$ (the Froissart-type bound this reflects is
the standard statement that a fixed-order, fixed-external-state
cross-section cannot rise without limit as $\sqrt s\to\infty$). A term
that grows like a positive power of $s$ at fixed, small $t_1,t_2$ — i.e.
fixed, small photon virtuality — is describing the near-real-photon corner
of phase space, where $e^-e^+\to e^-e^+\pi^0$ is dominated by
$\gamma\gamma$-fusion off nearly-on-shell photons radiated from the beams.
In that regime the process factorizes (equivalent-photon approximation)
into a photon flux (calculable, energy-dependent, but from ordinary QED
radiation — well understood and not the issue here) times the *real*
sub-process cross-section $\sigma(\gamma\gamma\to\pi^0)$, which is a fixed
resonance-formation cross-section that does **not** grow with the parent
beam energy. If an individual piece of the calculated amplitude survived,
uncancelled, growing as a power of $s$ in exactly this kinematic corner,
it would contradict that decomposition and imply an unbounded rate for a
process whose real physical sub-reaction has no such energy dependence —
a matter-of-principle violation, not a numerical nuisance. This is why the
passage calls it "problems with physics," not "problems with computing."

### Part 2: does the naive $s^2$ growth actually survive in $|\mathcal{M}|^2$?

**First check — derive the claimed $O(s^2M_\pi^4)$ numerator scaling
analytically**, using [`scripts/gram_s_scaling.py`](../scripts/gram_s_scaling.py).
At fixed small $t_1=-q_1^2,\,t_2=-q_2^2$ and $\sqrt s\to\infty$, the beam
momenta $p_1,p_2$ grow like $\sqrt s$ while the (near-real) photon momenta
$q_1,q_2$ stay $O(1)$ (their transverse component is fixed by $t_1,t_2$;
their longitudinal/energy components are suppressed by an extra power of
$1/\sqrt s$). Building the Gram (dot-product) determinant of
$(p_1,q_1,p_2,q_2)$ symbolically in this limit gives

$$\det\begin{pmatrix}p_1\!\cdot\!p_1 & p_1\!\cdot\!q_1 & p_1\!\cdot\!p_2 & p_1\!\cdot\!q_2\\ \cdots\end{pmatrix}
= -\tfrac14\,q_{T1}^2\,q_{T2}^2\,(\sqrt s)^4 = -\tfrac14\,t_1t_2\,s^2$$

i.e. $\varepsilon^{p_1q_1p_2q_2}\varepsilon_{p_1q_1p_2q_2}\sim s^2\,t_1t_2$
at fixed $t_1,t_2$ — confirmed by the script's `assert` on the leading
power of $\sqrt s$. Combined with `pi0.c`'s `part1 =
-64*levi.gram/(t_1t_2)^2`, this gives $\text{part1}\sim s^2/(t_1t_2)$: it
genuinely grows like $s^2$ at fixed $t_1,t_2$, exactly matching the
passage's $O(s^2M_\pi^4)$ claim for "the numerator... written in terms of
dotproducts."

**Second check — does this growth cancel in the full sum
`part1+part2+part3+part4`?** This is the part that needed direct numerical
verification against the real code, not assumption. Using
[`kinc2-driver/scan_pi0_scaling.c`](../kinc2-driver/scan_pi0_scaling.c) — an
instrumented driver built directly against the vendored `pickin.c`,
`orient.c`, `pi0.c` (see [`../kinc2/`](../kinc2/)), *not* a
reimplementation — $t_1,t_2$ are held fixed at $-0.001\,\text{GeV}^2$
(bisecting on the random numbers `pickin.c`'s `mapt()` uses to select
them, since their allowed range depends on $s$) while $\sqrt s$ is scanned
from $0.5$ to $10^4\,\text{GeV}$:

```
     sqrt(s)              s             t1             t2          part1          part2          part3          part4            sum
         0.5           0.25         -0.001         -0.001         511821        13317.9        15677.8        405.051         541222
           1              1         -0.001         -0.001    1.17694e+07        63430.6        75791.9        405.051     1.1909e+07
           3              9         -0.001         -0.001    1.17591e+09         631146         761754        405.051    1.17731e+09
          10            100         -0.001         -0.001    1.55476e+11    7.24568e+06    8.77605e+06        405.051    1.55492e+11
          30            900         -0.001         -0.001    1.28373e+13    6.58089e+07    7.97892e+07        405.051    1.28374e+13
        3000          9e+06         -0.001         -0.001    1.29588e+21    6.61046e+11    8.01879e+11        405.051    1.29588e+21
       1e+04          1e+08         -0.001         -0.001    1.59995e+23    7.34519e+12    8.91008e+12        405.051    1.59995e+23
```

`part1` dominates at every scanned energy and the **sum tracks `part1`
almost exactly** — from $\sqrt s=0.5$ to $\sqrt s=10^4$, $s$ grows by a
factor $4\times10^8$, and the sum grows by a factor
$\sim3\times10^{17}\approx(4\times10^8)^2$: genuine $s^2$ growth, with
essentially **no cancellation** between `part1..part4`. `part4` (the
$\lambda$ term, $s$-independent at fixed $t_1,t_2,M_\pi$) stays pinned at
$405.05$ throughout, confirming it plays no role in taming the growth
either.

**This directly matches the paper's own statement, previously
mis-paraphrased in this repo's notes**: Vermaseren 1983, p.353, states
explicitly *"in its Levi-Civita form all 4 terms are positive so no
cancellations occur"* — i.e. the paper itself says there is **no
cancellation between `part1..part4`**. The numerical run above confirms
this is exactly what the real code does: no cancellation, and the sum
genuinely grows like $s^2$ at fixed small $t_1,t_2$.

### So where does the "decent," non-growing behavior come from?

**It doesn't come from evaluating `pi0.c` as written and summing its four
terms — that sum is $O(s^2/(t_1t_2))$, not $O(1/(t_1t_2))$, confirmed
above.** The passage's "only $1/(t_1t_2)$ is decent" / "the numerator must
become very small" statement is not a description of `pi0.c`'s output. It
is a **statement about integrated cross-sections**, not the pointwise
squared matrix element: the actual differential cross-section carries
additional phase-space-measure suppression (from the Jacobian
$dj=1/\sqrt{-\Delta_4}$ that `pickin.c` computes and multiplies in — see
[Topic 1](topic1-peaks-and-O2n.md)'s discussion of $\Delta_4$) and, at
fixed physical beam energy, $t_1,t_2$ cannot be taken to zero at fixed $s$
independently of the rest of phase space — the accessible range
$[t_{1\min}(s),t_{1\max}(s)]$ itself depends on $s$, coupling the "$t_1,t_2$
small" and "$s$ large" limits through the phase-space boundary rather than
letting them vary independently the way the scan above (deliberately, to
isolate the matrix-element behavior) did. Disentangling exactly how the
Jacobian and the phase-space boundary combine with `pi0.c`'s $s^2$-growing
numerator to produce a finite, unitarity-respecting total cross section is
**not yet done here** — flagged as an open thread below rather than
asserted.

### Correction to Topic 4

[Topic 4](topic4-pseudoscalar-coupling.md) states *"every term in the
squared matrix element... is a genuine piece of the Gram-determinant
family, with nothing to cancel"* and *"consistent with the paper's claim
that all four terms carry the same sign and there is nothing to cancel."*
That statement about `part1..part4` not cancelling **is correct** and is
now backed by the numerical scan above (previously it was inferred from
the paper's prose alone). What was **not** in Topic 4, and needed this
topic to add: the fact that "no cancellation among `part1..part4`" is
entirely consistent with, and indeed requires, that the *raw* sum grows
like $s^2$ — the "decent $1/(t_1t_2)$" behavior the lecture note promises
is a property of the full cross-section (matrix element times phase-space
Jacobian times whatever additional $s,t_1,t_2$ structure survives
integration), not of `pi0.c`'s bare output taken alone.

### Open threads / not yet resolved

- The exact mechanism by which `pickin.c`'s Jacobian $dj$ and the
  $s$-dependent phase-space boundary combine with `pi0.c`'s $s^2$-growing
  numerator to produce a finite, non-power-growing physical cross-section
  has not been derived or checked numerically here — only the matrix
  element's own scaling has been isolated and confirmed.
- The "very bad cancellations between the various terms" that the passage
  warns about, if one naively contracted the Levi-Civita tensors into dot
  products, has not been checked directly (i.e. actually performing that
  dot-product expansion symbolically and confirming a large cancellation
  within it) — only inferred from `pickin.c`'s use of the stable factored
  form `levi.gram = (1-yy4^2)*dd/ap` (`pickin.c:188`), which structurally
  avoids ever forming that expansion.
