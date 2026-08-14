# Topic 4: why the pseudoscalar ($\pi^0$) coupling is the cleanest worked example

### The source passage

> "Let us study the reaction $e^-e^+\to e^-e^+\pi^0$. The $\pi^0$-particle is
> a pseudoscalar state and we will assume here for a moment that we do not
> have to worry about formfactors. We could also try a scalar state, which
> at the moment is more topical because of the Higgs particle, but the
> tricky point that we want to study shows itself more directly with a
> pseudo scalar state."
> — Vermaseren 1983, p.350

### The question

Everything established so far ([Topic 0](topic0-tree-vs-loop.md)–[Topic
3](topic3-subset-sum-proof.md)) is pure kinematics: the double-t-channel
$2\to3$ phase-space machinery (`pickin`/`orient`) that is identical
regardless of what particle sits at vertex $C$. So how is choosing a
*pseudoscalar* specifically related to "the tricky point" the paper wants
to study? This isn't a kinematics question at all — the answer lives in
**section 3 of the paper (the matrix element), not sections 1–2
(kinematics)**.

### What the pseudoscalar coupling actually forces

The $\pi^0\gamma\gamma$ coupling is $\mathcal{L}\sim\pi^0
F_{\mu\nu}\tilde F^{\mu\nu}=\pi^0\,\varepsilon^{\mu\nu\rho\sigma}
F_{\mu\nu}F_{\rho\sigma}$ — confirmed directly from the paper's own
formula (3.1):

$$|\mathcal{M}|^2 = g^2e^4\left|\frac{\bar u(p_3)\gamma_\mu u(p_1)\,
\varepsilon^{\mu\nu\rho\sigma}q_{1\rho}q_{2\sigma}\,
\bar v(p_2)\gamma_\nu v(p_5)}{q_1^2q_2^2}\right|^2$$

Parity forces this: there is **no other Lorentz structure** you can build
from two photon momenta $q_1,q_2$ and two photon polarization indices that
is a pseudoscalar. A scalar coupling doesn't have this constraint — the
paper says so directly (p.352): *"For a scalar one obtains similar
Levi-Civita tensors but there are also other terms."* Ordinary
parity-even dot-product structures are also available for a scalar, so
the amplitude is not forced to a single term.

### Why a single term matters: automatic, individual gauge invariance

The $\gamma^*\gamma^*\to\pi^0$ amplitude, $g\,\varepsilon^{e_1e_2q_1q_2}$,
is a **single term**, and — the paper states this explicitly (p.353,
directly under formula 3.2) — it is *"manifestly gauge invariant"*:
replacing either polarization vector $e_1\to q_1$ or $e_2\to q_2$ kills
the term identically, by the antisymmetry of $\varepsilon$ alone, with no
need for any other term to cancel against. Once this single term is
dressed with the electron/positron Dirac traces and squared, formula
(3.2) results:

$$|\mathcal{M}|^2 = g^2e^4[-64\varepsilon^{p_1q_1p_2q_2}
\varepsilon_{p_1q_1p_2q_2} - 16q_1^2\varepsilon^{p_2q_1q_2}
\varepsilon_{p_2q_1q_2\mu}-16q_2^2\varepsilon^{p_1q_1q_2\mu}
\varepsilon_{p_1q_1q_2\mu}-4q_1^2q_2^2\varepsilon^{q_1q_2\mu\nu}
\varepsilon_{q_1q_2\mu\nu}]/(q_1^2q_2^2)^2$$

— and the paper notes directly (p.353): *"If the Levi-Civita tensors are
contracted into 4-vector dot products this formula becomes numerically
unstable, but in its Levi-Civita form all 4 terms are positive so no
cancellations occur... the first term $\varepsilon^{p_1q_1p_2q_2}
\varepsilon_{p_1q_1p_2q_2}$ is the Gram determinant of the system and the
other 3 terms are minors of it, so all four terms have a meaning in the
kinematics."* Every term is a genuine, individually meaningful piece —
not an artifact of a gauge choice that must later cancel against another
term's artifact. This is directly visible in the actual code, `pi0.c`:

```c
double pi0(int par)
{
	double sum, tt = extra.t1*extra.t2;
	double la = (eemminput.mu*eemminput.mu-extra.t1-extra.t2)*0.5;
	double part1 = 0,part2 = 0,part3 = 0,part4 = 0;
	la = -(la*la-tt);
	part1 = -64*levi.gram/(tt*tt);
	part2 = -16*extra.t1*levi.dd2/(tt*tt);
	part3 = -16*extra.t2*levi.dd4/(tt*tt);
	part4 = -4*tt*la/(tt*tt);
	sum = part1+part2+part3+part4;
	return(sum);
}
```

Four terms, `part1..part4`, matching formula (3.2)'s four Levi-Civita
terms exactly, summed directly with **no subtraction anywhere in the
function** — consistent with the paper's claim that all four terms carry
the same sign and there is nothing to cancel.

### Why this makes $\pi^0$ the ideal *pedagogical* vehicle, not the hard case

The genuinely hard case — the one the whole paper (and the $O(2^n)$/
numerical-stability problem discussed in [Topic
1](topic1-peaks-and-O2n.md)) exists to solve — is
$\gamma^*\gamma^*\to\mu^+\mu^-$: **two separate Feynman diagrams**
(crossed muon-exchange, matching the box-diagram topology worked through
in [Topic 0](topic0-tree-vs-loop.md)'s Worked example 2), each of which
is gauge invariant only *jointly with the other*, not individually. Naive
4-vector-dot-product evaluation of either diagram alone contains large
terms that must cancel against the other diagram's terms to enforce
overall gauge invariance — precisely the "severe numerical cancellations"
problem stated in the paper's introduction (p.348: *"the very strong
gauge cancellations in the matrix element"*). Deriving that matrix
element (paper's eq. 3.3 onward: the identity introducing $K_1,K_2$
propagator denominators, the "super formulae," the four separate pieces
$M_{11},M_{12},M_{21},M_{22}$) needs substantially more machinery than
the $\pi^0$ case.

So: **the pseudoscalar choice is not about the reaction being harder or
requiring different kinematics — it's the simplest possible test case for
the Levi-Civita/Gram-determinant *method*** the paper needs for the
genuinely hard multi-diagram QED case. Parity forces $\gamma\gamma\to\pi^0$
into a single, automatically gauge-invariant term, so Vermaseren can
demonstrate "write everything in Levi-Civita form instead of dot
products, no cancellations occur" cleanly, on the simplest possible
example, before tackling the case where it's actually needed to solve a
real numerical problem.

See [Topic 5](topic5-s-scaling-pi0.md) for the follow-up question this
raises: "no cancellation between `part1..part4`" does not by itself mean
the sum is well-behaved in $s$ — checked directly, the sum grows like
$s^2$ at fixed small $t_1,t_2$, exactly as much as `part1` alone.

### Open threads / not yet resolved

- We have not read the paper's actual derivation of the muon-pair matrix
  element (section 3, eq. 3.3 onward, the "super formulae") in detail —
  only confirmed *why* it's harder than the $\pi^0$ case, not worked
  through *how* the method handles the two-diagram gauge-invariance
  problem concretely.
- The paper's remark that scalar-coupling's extra terms "do not
  contribute to the major numerical disasters" (p.352) is stated but not
  derived here — we have not checked why those particular extra terms
  happen to be numerically benign, only noted that they exist (unlike the
  pseudoscalar case, where no such extra terms exist at all).
