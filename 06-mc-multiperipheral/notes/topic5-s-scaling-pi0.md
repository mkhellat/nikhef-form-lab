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

   **(Corrected 2026-11: "$s$-independent" above is too strong. The
   integrated cross-section is *not* $s$-independent: without cuts it grows
   like $\ln^3 s$. What the $1/(t_1t_2)$ behaviour buys is growth that is
   *only logarithmic*, never a power of $s$. See
   [The integrated law, resolved](#the-integrated-law-resolved-2026-11)
   below.)**

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

**(Corrected 2026-11: the argument above is right about *powers* of $s$,
but this repo read it as "the integrated cross-section must not grow at
all" (see the "rather than flattening to a constant" expectation further
down). That is too strong. Two things were glossed over:**

- **The photon flux grows.** The paragraph itself says the flux is
  "energy-dependent". Each beam radiates an equivalent-photon spectrum
  whose size grows like $\ln s$ (the photon can be emitted over a wider and
  wider range of virtualities and energy fractions). So
  $\sigma(e^+e^-\to e^+e^-\pi^0)$ = flux $\otimes$ flux $\otimes$ (fixed
  $\sigma(\gamma\gamma\to\pi^0)$) grows like a polynomial in $\ln s$, and
  with no cuts at all it grows like $\ln^3 s$. That is derived and checked
  against the real code below. A non-growing $\sigma(\gamma\gamma\to\pi^0)$
  does *not* mean a non-growing $e^+e^-$ cross-section.
- **"Cannot rise without limit" is not true for this fixed-order rate.**
  The uncut tree-level rate computed here *does* rise without limit, but
  only logarithmically. The Froissart-type remark should be read loosely:
  what unitarity-style reasoning rules out is growth like a *power* of $s$
  that the factorised picture cannot explain. Log growth that comes from
  the fluxes is expected. It is what the lecture means by a "decent"
  cross-section (part5.pdf p.1: "Only a behaviour like $1/(t_1t_2)$ could
  give a decent crosssection").

**The original paragraph is left intact above so the reasoning trail stays
visible.)**

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

**(Resolved 2026-11: the Jacobian does the job, and the boundary does not.
`pickin.c`'s Jacobian carries $1/s$ and $1/s_2$ factors. Together with the
$1/(2s)$ flux they cancel the $s^2$ of `part1` exactly, which leaves a
density that is flat in $(\ln|t_1|,\ln|t_2|,\ln s_2)$. The integral is
then a phase-space *volume* in log variables, and it grows like $\ln^3 s$,
not like a power. See
[The integrated law, resolved](#the-integrated-law-resolved-2026-11).)**

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

**(Note 2026-11: the "whatever additional structure" is now pinned down.
The Jacobian $\times$ flux removes the $s^2$, leaving $(dt_1/t_1)(dt_2/t_2)
(ds_2/s_2)$. "Decent" means *logarithmic* growth, not a constant. See
[The integrated law, resolved](#the-integrated-law-resolved-2026-11).)**

### The Jacobian/boundary mechanism, partially resolved

The obvious first guess — that the accessible $t_1,t_2$ range shrinks as
$s\to\infty$, squeezing out the growing region — is **wrong** in this
process. Checked directly (massless $e^-/e^+$ beams, $m_1=m_2=m_3=0$):
`pickin.c`'s own $t_{1\max}$ formula (`pickin.c:90`) simplifies exactly to
$t_{1\max}=\sigma_1-s\to -s$, and $t_{1\min}=0$ identically in this mass
limit (both confirmed against the running code, not just algebra). So the
kinematically allowed range is $t_1\in[-s,0]$ — it **grows** with $s$, it
does not shrink; the near-real-photon corner ($t_1\to0^-$) remains
accessible at every $s$, at all energies, exactly where `part1`'s $s^2$
growth lives.

So the resolution has to come from the *integral*, not the boundary
alone: the true cross section is $\sigma\propto\int dj\,|\mathcal M|^2$
over the full $s$-dependent phase space (`kinc2/eee.c`'s driver,
`eepi.c`, computes exactly `overallconstant * dj * pi0(0)` per phase-space
point — read directly from the un-vendored reference archive, not
guessed), where `eee.c` shows `overallconstant` carries an explicit
$1/(2s)$ flux factor (the standard 2-body-flux normalization) on top of
$s$-independent couplings.
[`kinc2-driver/mc_xsec_s_scan.c`](../kinc2-driver/mc_xsec_s_scan.c) Monte
Carlo-integrates $dj\times\text{pi0}(0)$ over the full unit hypercube (no
importance sampling — a crude, unbiased estimate, not a precision VEGAS
run) at several $\sqrt s$, then applies that same $1/(2s)$ factor:

```
   sqrt(s)          raw_avg     raw_avg/(2s)
         1     1.093544e+04     5.467718e+03
         3     2.671909e+05     1.484394e+04
        10     6.345207e+06     3.172604e+04
        30     9.747381e+07     5.415212e+04
       100     1.752242e+09     8.761212e+04
```

**(Corrected 2026-11: "no importance sampling, a crude estimate" was
wrong. Uniform points on the unit cube are *not* uniform in $t_1,t_2,s_2$.
`pickin()` turns each random number into a kinematic variable with a
logarithmic map:**

- **$t_1$ and $t_2$ use `mapt(tmin,tmax,&dt,x,1)` (`kinc2/pickin.c:93`,
  `:142`). Type 1 is commented "Assumes dt/t" (`kinc2/mapt.c:13-20`): it
  sets $t=t_{\min}(t_{\max}/t_{\min})^x$ with weight
  $dt=-t\ln(t_{\max}/t_{\min})$. So $x$ is uniform in $\ln|t|$.**
- **$s_2$ uses `mapw(sig1,smax,&ds2,x,1)` (`kinc2/pickin.c:79`, the
  `option == 0` branch that `orient()` takes, `kinc2/orient.c:31`). Type 1
  is "Assumes dw2/w2" (`kinc2/mapw.c:13-20`), so $x$ is uniform in $\ln s_2$.**

**That is importance sampling, and for this integrand it is close to ideal.
The real density is flat in exactly these log variables (derived below),
so $f(r)$ is nearly constant on the cube. That is why plain random points
already work well. The merged VEGAS driver
[`kinc2-driver/xsec_pi0_vegas.c`](../kinc2-driver/xsec_pi0_vegas.c)
(header, lines 6-13) reproduces the numbers above to within 0.01-0.1%.
Re-run 2026-11, `raw_avg/(2s)` above vs VEGAS (200k calls x 10 iterations):
$\sqrt s=1$: 5467.7 vs 5464.8; 10: 31726.0 vs 31695.8; 100: 87612.1 vs
87654.1. Its `--check` gate G2 pins this agreement at 0.3%. The table is
kept as the original record.)**

The raw integral (before the flux factor) grows with a *decreasing*
effective power of $s$ as $s$ increases (fit between consecutive points:
$s^{1.45}\to s^{1.31}\to s^{1.24}\to s^{1.20}$), consistent with
integrated growth closer to $s\log s$ than to `part1`'s bare $s^2$ — the
phase-space integration measure and the boundary do suppress the growth
substantially, as the lecture claims. But after dividing by the $1/(2s)$
flux factor, the result **still grows** — roughly $\times16$ from
$\sqrt s=1$ to $\sqrt s=100$ — rather than flattening to a constant. This
is a genuine, only partially resolved result, not a derivation error to
paper over: the *unrestricted*, full-$4\pi$ acceptance integral computed
here plausibly diverges logarithmically as $\sqrt s\to\infty$ (a known
feature of equivalent-photon-type cross sections when integrated over
photon virtuality all the way down to the exactly-real-photon point,
regularized physically by the electron mass rather than by any dynamics
in $\pi^0\gamma\gamma$ coupling) — but that asymptotic log-law has not
been derived analytically here, only observed numerically over a modest
range of $s$ with a statistically noisy MC estimate. A real experiment's
cross section is finite because detector acceptance always imposes a
nonzero minimum scattering angle (equivalently, a nonzero minimum
$|t_1|,|t_2|$), which directly excludes the region `part1`'s $s^2$ growth
comes from; that acceptance-cut version of the cross section was not
computed here.

**(Resolved 2026-11, one sentence at a time. The paragraph above is kept
as written.)**

- **"plausibly diverges logarithmically ... has not been derived
  analytically": now derived.** After the flux factor, the uncut integral
  grows like $\ln^3 s$. Before it, the raw integral grows like
  $s\ln^3 s$, not the "$s\log s$" guessed above. The closed-form
  leading-log law is given in the next section. The
  "decreasing effective power" $s^{1.45}\to s^{1.20}$ and the "roughly
  $\times16$" are both explained by that law plus $O(1/\Lambda)$
  corrections.
- **"statistically noisy MC estimate": wrong** (see the correction above).
  The numbers are good to about 0.1%. The scan now runs up to
  $\sqrt s=10^6\,$GeV with VEGAS.
- **"a nonzero minimum scattering angle ... excludes the region": half
  right.** A minimum *angle* does give a flat cross-section, about 1300 in
  driver units from 30 GeV to 100 TeV. A fixed minimum $|t|$ does *not*:
  the cross-section still grows. The reason is that with a pointlike
  $\pi^0$ vertex, the growth does not come only from the $t\to0$ corner. It
  also comes from the tail out to $|t|\sim s$. See
  [Acceptance cuts](#acceptance-cuts-what-makes-the-cross-section-decent).

### The integrated law, resolved (2026-11)

**Executable evidence:** `make -C 06-mc-multiperipheral/kinc2-driver
check-xsec-pi0`. It runs
[`kinc2-driver/xsec_pi0_vegas.c`](../kinc2-driver/xsec_pi0_vegas.c), a
VEGAS integration of the same integrand
$f(r)=dj(r)\,|\mathcal M|^2_{\pi^0}/(2s)$ built from the vendored
`orient()`/`pickin()`/`pi0()`. Its header (lines 55-77) states the LL laws
and its check gates (lines 79-88). Output is in "raw/2s" units. Only the
*shape* in $s$ is meaningful, because the absolute normalisation is still
open (see the last subsection). Re-run 2026-11: PASS in about 11 s.

#### Step 1: why the density is flat in log variables

Three ingredients, each read off the code:

1. **The matrix element.** Vermaseren treats the $\pi^0$ "like a point
   particle, i.e. without form-factors" (Vermaseren 1983, p.352). So
   nothing in $|\mathcal M|^2$ switches off at large $|t|$. Part 2 above
   showed $|\mathcal M|^2\approx\text{part1}=-64\,\text{gram}/(t_1t_2)^2$
   with $\text{gram}\simeq-\tfrac14 t_1t_2s^2$. So
   $|\mathcal M|^2\sim s^2/(t_1t_2)$: this is the lecture's
   "$1/(t_1t_2)$" behaviour, with the $s^2$ that Part 2 found.
2. **The Jacobian.** `pickin.c:187` returns
   `dj = ds2*dt1*dt2*pi*pi/(8*pickzz.sl1*sqrt(-ap))`. Here
   `sl1` $=\lambda^{1/2}(s,m_e^2,m_e^2)\approx s$ (`pickin.c:71,75`), and
   `ap` $=-\tfrac14\lambda(s_2,t_1,m_e^2)$ (`pickin.c:170`), so
   $\sqrt{-ap}\approx s_2/2$ once $s_2\gg|t_1|$. Here $s_2$ is the
   invariant mass² of the $\pi^0$+positron system. So
   $dj\propto dt_1\,dt_2\,ds_2/(s\,s_2)$.
3. **The flux** $1/(2s)$, from `eee.c`'s `overallconstant`
   (`kinc2/eee.c:27-28` in the un-vendored archive, see
   [`../kinc2/README.md`](../kinc2/README.md)).

Multiply them:

$$d\sigma\;\propto\;\frac{1}{2s}\cdot\frac{dt_1\,dt_2\,ds_2}{s\,s_2}\cdot
\frac{s^2}{t_1t_2}\;=\;\frac{dt_1}{t_1}\,\frac{dt_2}{t_2}\,\frac{ds_2}{s_2}.$$

**The $s^2$ is gone.** The density per unit $\ln|t_1|$, $\ln|t_2|$, $\ln
s_2$ is a constant that does not depend on $s$. (At fixed $t$'s, $\ln
s_2$ is essentially the $\pi^0$ rapidity.) So the cross-section is
(constant) $\times$ (volume $V$ of the allowed region in these log
variables). It grows only because the log ranges widen as $s$ grows.
Those are exactly the variables `pickin` already samples uniformly (see
the correction above).

**Why $1/(t_1t_2)^2$ would have meant powers of $s$** (this makes the
lecture sentence concrete; it is a derivation, not a run). The smallest
allowed $|t_i|$ is about $m_e^2x_i^2$, where $x_i$ is photon $i$'s energy
fraction (`pickin.c:91-92` gives $t_{1\min}$ exactly). Making the $\pi^0$
needs $x_1x_2s\approx M^2$, so typically $x_i\sim M/\sqrt s$ and
$|t_{i,\min}|\sim m_e^2M^2/s$. Compare the two integrals:

- $\int dt/t$ gives $\ln(1/|t_{\min}|)$, which grows like $\ln s$.
  Harmless.
- $\int dt/t^2$ gives $1/|t_{\min}|\sim s/(m_e^2M^2)$. That is a power of
  $s$ per photon, and $s^2$ for both photons.

The numerator (the Gram determinant) vanishing like $t_1t_2$ is what turns
the second case into the first.

#### Step 2: the volume, in closed form

Notation: $M=m_{\pi^0}$, $A=\ln(M^2/m_e^2)\approx11.15$,
$\Lambda=\ln(s/M^2)$. Use $u_i=-\ln x_i$. Then $u_1+u_2=\Lambda$ (the
$\pi^0$ mass condition), and photon $i$'s $\ln|t_i|$ runs from
$\ln(m_e^2x_i^2)$ up to the region's upper edge.

- **$|t_1|,|t_2|<M^2$.** Each photon's log range is $A+2u_i$. So
  $V=\int_0^\Lambda du\,(A+2u)(A+2\Lambda-2u)=A^2\Lambda+2A\Lambda^2+
  \tfrac23\Lambda^3$. This is checked symbolically by sympy, in the
  scratch script `llvol.py` described under "Status of this evidence"
  below.
- **Uncut.** Above $|t|=M^2$ the pointlike vertex keeps contributing, out
  to $|t|\sim s$. The volume of the extra polytope gives
  $V=A^2\Lambda+3A\Lambda^2+\tfrac53\Lambda^3$, which grows like
  $\ln^3s$. This closed form matches a direct numerical polytope
  integration (`llvol.py`) to $10^{-5}$ for $\sqrt s=1$ to $10^4\,$GeV
  (re-run 2026-11).

These are the formulas in the driver header (`xsec_pi0_vegas.c:64-65`).
The method is the equivalent-photon/Weizsäcker-Williams picture. The
header cites Budnev et al., Phys. Rep. 15 (1975) 181 for it. That paper is
not in `__src__/_ref_/` and was not checked here.

**Effective power.** The 08-2026 power fits used `raw_avg`, which is
measured *before* the $1/(2s)$ flux, so `raw_avg` $\propto s\cdot\sigma$.
Its local power is
$p=1+d\ln\sigma/d\ln s=1+d\ln V/d\Lambda$. With $V\sim\Lambda^3$ this
gives $p\to1+3/\Lambda$. The power drifts down toward 1 forever and never
settles on a constant. That is the recorded "$s^{1.45}\to s^{1.2}$".

#### Step 3: MC vs LL

The rows come from `compare_ll.py`, run on VEGAS data (scratch, see
below). The $\sqrt s=1,10,100,1000$ rows were re-confirmed with the merged
driver in 2026-11. $p$ is computed between consecutive rows. MC/LL uses
the LL constant $\pi^2$ (raw/2s per unit $V$), which assumes the factor ½
discussed below. Its *level* is therefore provisional. Its *trend* is not.

```
  sqrt(s)/GeV  Lambda   sigma_raw(MC)   MC/LL   p_MC    p_LL
        1       4.00    5.4649e+03     0.485     -       -
        3       6.20    1.4845e+04     0.612   1.455   1.349
       10       8.61    3.1701e+04     0.696   1.315   1.262
       30      10.81    5.4114e+04     0.745   1.243   1.212
      100      13.22    8.7642e+04     0.783   1.200   1.180
     1000      17.82    1.8256e+05     0.830   1.149   1.138
     1e6       31.64    7.9590e+05     0.894   1.090   1.086
```

How to read it:

- $p_{\rm MC}$ follows $p_{\rm LL}$ downward at every point.
- MC/LL rises steadily toward 1. That is the signature of $O(1/\Lambda)$
  corrections, which the leading-log counting drops. (They come from
  $s_2\gg|t|$ not being exact, the $(1-x)$ factors, and the edges of the
  region.)
- The "$\times16$ from 1 to 100 GeV" is MC 16.04, while LL alone gives
  9.9. At low $\Lambda$ the subleading terms are large and increase the
  growth. Asymptotically the two agree better: $\sigma(10^6)/\sigma(10^3)$
  is 4.36 (MC) vs 4.05 (LL).

#### Acceptance cuts: what makes the cross-section "decent"

- **Fixed $|t_1|,|t_2|\ge t_0$: still grows.** The LL volume is
  $V=B^2\Lambda+B\Lambda^2+\Lambda^3/3$ with $B=\ln(M^2/t_0)$
  (`xsec_pi0_vegas.c:66`). Removing the $t\to0$ corner removes the $A$
  terms, but the pointlike tail up to $|t|\sim s$ still widens like
  $\Lambda$. VEGAS with $t_0=0.01\,$GeV² gives MC/LL $=0.36\to0.85$ for
  $\sqrt s=1\to10^5\,$GeV, again rising toward 1 (scratch data
  `build/t001.txt`).
- **Fixed electron angle $\theta\ge\theta_{\min}$: flat.** Dimensional
  argument: at fixed angle, the smallest allowed $|t|$ is about
  $s\theta^2/4\cdot(1-x)$, and the largest is about $s$. So the
  $\ln|t|$ window has an $s$-independent width $\ln(4/\theta^2)$. Once
  $M^2\ll s\theta^2/4$, i.e. $\sqrt s\gg2M/\theta$ (9 GeV for 30 mrad),
  no remaining log range depends on $s$ and $V$ goes to a constant. This
  is an argument, not a derivation of the constant. The naive LL
  asymptote $\pi^2\ln^3(4/\theta^2)/3\approx1950$ does not match the MC's
  1300. **Measured:** gate G5 (`xsec_pi0_vegas.c:80`), 30 mrad, re-run
  2026-11: $\sigma_{\rm raw}=1298.3,\ 1300.8,\ 1301.2$ at
  $\sqrt s=30,100,1000\,$GeV. A separate re-run gave 1299.9, 1300.7,
  1297.9 ($\pm0.7$) at $10^3,10^4,10^5\,$GeV. **Flat at about 1300 from
  30 GeV to 100 TeV.** This is the concrete "decent cross section" check
  the old open thread asked for.
- **With a physical form factor** (not modelled in `pi0.c`, which is
  pointlike), large $|t|$ would be suppressed. Under a fixed angle cut,
  the minimum $|t|$, about $s\theta^2/4$, then runs into the form-factor
  scale, and $\sigma$ would *fall* with $s$. This is an argument only, not
  computed.

#### Where the real "very bad cancellations" are

They are not between `part1..part4` (Part 2 above; Vermaseren p.353). They
are *inside* the Gram determinant once it is expanded into dot products.
The lecture's "terms that are very big ($O(s^2M_\pi^4)$)" are individual
products such as $(p_1{\cdot}p_2)^2(q_1{\cdot}q_2)^2\approx(s/2)^2(M^2/2)^2$.
The sum of all such products is $\simeq-\tfrac14t_1t_2s^2$. So the
expansion cancels by a factor of about $M^4/(t_1t_2)$, which does not
depend on $s$ but becomes huge as $t\to0$, exactly where the cross-section
lives. (Throwaway mpmath check, 2026-11, not committed: at
$t_1=t_2=-10^{-4}\,$GeV², the largest Leibniz term is $\approx s^2M^4/16$
and exceeds the determinant by $1.2\times10^4$ at $\sqrt s=10,100,1000$.
For a general azimuth $\phi$ between the two photons, the determinant is
$-\tfrac14t_1t_2s^2\sin^2\phi$. Part 2's script uses $\phi=90°$.)

`pickin.c:188` never forms that expansion. It builds the Gram determinant
in factorised form, `levi.gram = (1-yy4)*(1+yy4)*dd/ap`. Here
`yy4 = cos(pi*rannums[3])` is the cosine of the azimuthal variable that
`pickin` samples (`pickin.c:177`), so $1-y^2$ is a $\sin^2$. `dd` is a
product of boundary-root factors $(s_2-s_{2\max})(s_2-s_{2\pm})\cdots$
(`pickin.c:171-178`). Each factor is computed directly, so no large terms
are ever subtracted. [Topic 4](topic4-pseudoscalar-coupling.md) (section
on the Gram minors, "29 times too large")
already demonstrates the same mechanism numerically, for the `part3` minor
with [`check_gram_minors.c`](../kinc2-driver/check_gram_minors.c): the
naive minor comes out 29$\times$ too large and with the wrong sign at
$t_1\sim10^{-11}$. That is not repeated here.

#### Status of this evidence

- **In the repo:** the driver and its `check-xsec-pi0` gates G2, G4, G5.
- **Scratch only, not committed:** the LL volume derivation (`llvol.py`,
  sympy + mpmath), the MC/LL table (`compare_ll.py`), the EPA comparison
  (`epa.py`) and the VEGAS scan data (`build/*.txt`). These came from an
  earlier analysis pass and were re-run 2026-11. If they are worth
  keeping, they belong in `scripts/` (separate commit).

#### Still open: absolute normalisation

Two things are still unresolved.

1. **The absolute $\pi^0$ normalisation.** First put back the
   $(2\pi)^{-5}$ phase-space factor that `eee.c`'s `overallconstant`
   lacks (`kinc2/eee.c:27-28`). Then compare kinc2 with a numerical
   equivalent-photon (EPA) integral, restricted to
   $|t_1|,|t_2|<M^2$, the region where the EPA with $Q^2_{\max}=M^2$ is
   the right comparison. The ratio kinc2/EPA is 0.475, 0.491, 0.495 and
   0.497 at $\sqrt s=1,10,100,1000\,$GeV (`epa.py`, re-run 2026-11). That
   is consistent with exactly ½ asymptotically, approached as
   $O(1/\Lambda)$. It is not shown to be exactly ½ at finite $s$. Where
   the factor 2 belongs is unresolved. It could be `eee.c`'s
   $g^2=7.8\times10^{-9}\cdot32\pi/M^3$ (`eee.c:21`), or `pi0.c`'s
   spin-sum/averaging convention. The uncut comparison against an EPA with
   $Q^2_{\max}\sim s$ gives 0.36 to 0.28 and is *not* a clean test,
   because the transverse-photon EPA formula does not model the pointlike
   $|t|\sim s$ tail. A missing factor 2 is a constant, so it cannot
   change any of the scaling conclusions above. Until it is settled, the
   driver prints no picobarns (`xsec_pi0_vegas.c:36-42`).
2. **The $O(1/\Lambda)$ coefficients** behind the MC/LL ratio rising
   toward 1, and the constant of the angle-cut plateau, have not been
   derived. Both are known only numerically.

### Open threads / not yet resolved

- The exact asymptotic power-or-log law of the fully-inclusive (uncut)
  integrated cross section as $\sqrt s\to\infty$ has not been derived
  analytically — only estimated numerically over $\sqrt s=1$–$100\,$GeV
  via a crude, unimportance-sampled Monte Carlo integral, which is too
  noisy and too narrow a range to distinguish a slow power law from a
  logarithm with confidence.
  **(Resolved 2026-11: the law is a cubic polynomial in $\Lambda=\ln(s/M^2)$,
  $V=A^2\Lambda+3A\Lambda^2+\tfrac53\Lambda^3$ at leading log, so it grows
  like $\ln^3s$ and is not a power. It is confirmed by VEGAS up to
  $\sqrt s=10^6\,$GeV, with MC/LL rising $0.49\to0.89$ toward 1. The
  premise was also wrong: the MC was importance-sampled in $dt/t$
  (`pickin.c:93,142`) and $ds_2/s_2$ (`pickin.c:79`), and it agrees with
  VEGAS to about 0.1%. See
  [The integrated law, resolved](#the-integrated-law-resolved-2026-11).)**
- The finite, $t_1,t_2$-acceptance-cut cross section (the physically
  realistic quantity any real detector measures) was not computed —
  doing so and checking that IT is flat/decreasing in $s$ (unlike the
  uncut integral above) would be the more direct confirmation of the
  lecture's "decent cross section" claim.
  **(Resolved 2026-11: a fixed electron angle $\theta\ge30\,$mrad gives
  $\sigma_{\rm raw}\approx1300$, flat from 30 GeV to 100 TeV, checked by
  gate G5 of `make -C 06-mc-multiperipheral/kinc2-driver check-xsec-pi0`.
  A fixed $|t|\ge t_0$ cut does *not* give a flat result. It still grows,
  because the pointlike tail extends to $|t|\sim s$. See
  [Acceptance cuts](#acceptance-cuts-what-makes-the-cross-section-decent).)**
- **Still open (2026-11):** the absolute $\pi^0$ normalisation. With
  $(2\pi)^{-5}$ restored (missing from `eee.c:27-28`), kinc2 comes out at
  ½ of the EPA prediction (0.475 rising to 0.497 over 1 GeV to 1 TeV, with
  $|t|<M^2$). It is not settled whether the factor 2 belongs in `eee.c`'s
  `g2` or in `pi0.c`'s spin sum. It does not depend on $s$, so the
  scaling results are unaffected.
- **Still open (2026-11):** the $O(1/\Lambda)$ corrections to the LL law
  (the MC/LL ratio is 0.61 at 3 GeV and 0.83 at 1 TeV) and the value of
  the angle-cut plateau (the naive LL estimate of about 1950 vs the
  measured 1300). Neither has been derived.
- **Still open (2026-11):** the LL scripts (`llvol.py`, `compare_ll.py`,
  `epa.py`) exist only as scratch files. Committing them under
  `scripts/` would make the closed-form volumes and the EPA ratio
  reproducible from the repo, as the driver already makes the scaling
  reproducible.
