# Topic 2: "the number of non-trivial integration variables is at most $3n-4$"

### What "non-trivial integration variable" means

Phase space for an $n$-particle final state is, before any reduction,

$$d\Phi_n = \prod_{i} \frac{d^3p_i}{(2\pi)^3\,2E_i}\;(2\pi)^4\,
\delta^{(4)}\!\Big(P-\sum_i p_i\Big)$$

Each final-state particle contributes 3 real numbers (its 3-momentum —
energy is fixed by the mass-shell condition $E_i=\sqrt{p_i^2+m_i^2}$, so it
is *not* an independent degree of freedom). That gives $3n$ numbers before
any constraint. The $\delta^{(4)}$ removes exactly 4 of them (fixes total
energy and 3-momentum to match the initial state). Nothing else is
subtracted in the formula itself:

$$3n - 4$$

**"Non-trivial" here means: an independent number you actually integrate
(sample) over**, as opposed to (a) a component fixed outright by
conservation, or (b) a genuinely flat/uncorrelated direction that
contributes only an overall constant to the integral (e.g. $2\pi$) without
adding any structure worth stratifying/mapping.

**Now checked directly against the primary reference, not just derived
independently.** Byckling & Kajantie, *Particle Kinematics* (1973),
Table III.1 (p.49) — the standard reference this whole counting scheme
traces back to — gives *four* distinct variable counts for a $2\to n$
process, not one, and it matters which one "$3n-4$" refers to:

| Byckling–Kajantie row | count ($2\to n$) | example ($2\to2$) |
|---|---|---|
| All Variables | $3n-3$ | $s,t,\phi$ |
| Essential Variables | $3n-4$ | $s,t$ |
| **Final State Variables** | $3n-4$ | $t,\phi$ |
| Essential Final State Variables | $3n-5$ | $t$ |

Two different rows both give $3n-4$, for different reasons: "Essential
Variables" ($s,t$) counts the whole process including $s$ (the total
CM energy², fixed by the beam energy — not something a Monte Carlo run
integrates over) while already dropping $\phi$ as trivial; "**Final
State Variables**" ($t,\phi$) counts only the outgoing particles'
freedom (correctly excluding $s$, since $\sqrt{s}$ is a fixed
experimental input here, matching `eepi.c`'s `eemminput.sq`) while
*including* $\phi$. Vermaseren's sentence, and this document's Topic 2,
are about final-state phase-space integration for a fixed beam
energy — so the relevant row is **"Final State Variables," $3n-4$,
which does include the trivial $\phi$** — matching what this document
independently derived from the code below. The further reduction this
document found (`pickin`/`orient`/`eepi` sample only 4 of the formal 5
slots for $n=3$, dropping $\phi$ as computed-not-sampled) is exactly
Byckling–Kajantie's own next row, "**Essential Final State Variables**,"
$3n-5$ — for $n=3$: $3n-5=4$, matching exactly.

**Check against the familiar 2-body case ($n=2$):** "Final State
Variables" $=3n-4=2$: $(\theta,\phi)$, or equivalently $(t,\phi)$ per the
table's own example, since $\theta$ and $t=(p_1-p_3)^2$ carry the same
information. $\theta$/$t$ is physical. $\phi$ (azimuth around the beam
axis) is "trivial" in the sense above: for unpolarized beams the physics
doesn't depend on it, it integrates to a flat factor of $2\pi$ — but it
is still formally one of the $3n-4$ "Final State Variables," not a
*third* subtraction on top of $3n-4$; Byckling–Kajantie's own "Essential
Final State Variables" row performs that subtraction separately, giving
$3n-5=1$ (just $t$) for $n=2$. This matters: **the naive expectation
that you should separately subtract 3 more for "overall rotational
freedom of the final state" is wrong** — only 1 direction (azimuth about
the fixed beam axis) is a symmetry for a fixed-beam collision; the polar
orientation relative to the beam is physical. So $3n-4$ ("Final State
Variables") already has exactly the right amount of "triviality" baked
in via the conservation-law subtraction alone, and the one further
trivial direction ($\phi$) is removed by a *separate*, well-defined next
step ($3n-4\to3n-5$) rather than needing 3 more ad hoc subtractions.

### Checking this against our actual code: `pickin.c`/`orient.c`/`eepi.c` ($n=3$)

$3n-4 = 3(3)-4 = 5$ for our diagram. Every claim below is anchored to
specific lines, re-extracted and read directly from `kinc2.tar.gz`
(`_course_/-0-uam19/kinc2.tar.gz`) on 2026-08-08 — quoted, not paraphrased,
so this can be independently re-checked against the archive.

**Claim: only 4 random numbers are consumed by the $n=3$ pipeline, not 5.**

`eepi.c:14-20` is the entry point VEGAS would call for this reaction:

```c
double eepi(double *rannums)
{
    ...
    dj = orient(eemminput.sq*eemminput.sq,
        eemminput.me,eemminput.me,eemminput.me,eemminput.mu,eemminput.me,
        rannums,0);
```

`orient.c:31` immediately forwards the same `rannums` pointer into `pickin`
unchanged:

```c
if ( ( dj = pickin(s,m1,m2,m3,m4,m5,rannums,options) ) == 0 ) return(0);
```

Inside `pickin.c`, a full-text search for `rannums[` turns up exactly these
four indices and no others (`pickin.c:79,93,126,127,142,167,168,177`):

```c
79:   extra.s2 = mapw(sig1,smax,&ds2,rannums[2],1);
93:   extra.t1 = mapt(t1min,t1max,&dt1,rannums[0],1);
126:  if ( option < -1 ) extra.s2 = mapw(sig2,s2max,&ds2,rannums[2],1);
127:  if ( option == -1 ) extra.s2 = mapla(sig2,s2max,extra.t1,pickzz.w2,&ds2,rannums[2]);
142:  extra.t2 = mapt(t2min,t2max,&dt2,rannums[1],1);
167:  if ( option > 1 ) extra.s2 = mapw(s2min,s2max,&ds2,rannums[2],1);
168:  if ( option == 1 ) extra.s2 = mapla(s2min,s2max,extra.t1,pickzz.w2,&ds2,rannums[2]);
177:  yy4 = cos(pi*rannums[3]);
```

So the indices actually touched are `{0,1,2,3}` — never `rannums[4]`. This
is independently corroborated by `eepi.c`'s own diagnostic/histogram block
(`eepi.c:29-39`), which logs exactly these four back out and nothing more:

```c
if ( withplots ) {
    lplot.xvalue[0] = log10(-extra.t1);
    lplot.xvalue[1] = log10(-extra.t2);
    lplot.xvalue[2] = log10(variab.p3*variab.st3);
    lplot.xvalue[3] = log10(variab.p4*variab.st4);
    lplot.xvalue[4] = log10(variab.p5*variab.st5);
    lplot.xvalue[5] = rannums[0];
    lplot.xvalue[6] = rannums[1];
    lplot.xvalue[7] = rannums[2];
    lplot.xvalue[8] = rannums[3];
}
```

**Caveat on this claim:** neither `kinc.tar.gz` nor `kinc2.tar.gz` ships a
`main`-style driver that calls `vegas(eepi, ..., ndim, ...)` for this
specific reaction (unlike e.g. `kinc/main3.c:70`:
`vegas(feemm,(double)0.00001,2,60000,6,1)`, which shows `feemm` is called
with VEGAS dimension 6 — that pattern would let us read off `ndim`
directly if a `runpi`-equivalent main file existed). So "4 random numbers,
hence VEGAS dimension 4" is an *inference from usage inside `pickin`*, not
a directly observed `vegas(eepi,...)` call. It is a strong inference (a
function can't productively read an index it doesn't touch), but it should
be flagged as such rather than stated as if we'd seen the `vegas()` call
itself — we haven't, for this particular reaction.

**Claim: what the 4 random numbers become, exactly.**

Reading `pickin.c` top to bottom, in call order:

- `rannums[0]` → `pickin.c:93`, `extra.t1 = mapt(t1min,t1max,&dt1,rannums[0],1)` — the first propagator invariant $t_1$.
- `rannums[1]` → `pickin.c:142`, `extra.t2 = mapt(t2min,t2max,&dt2,rannums[1],1)` — the second propagator invariant $t_2$, computed *after* `t1` because `t2min`/`t2max` depend on `extra.t1` (see `d4 = pickzz.w4-extra.t1` at `pickin.c:96`, used in the `t2max`/`t2min` formulas at lines 139-141). This ordering dependency is itself worth remembering: the 4 variables are not sampled independently of each other's *ranges*, only independently as [0,1]-uniform inputs.
- `rannums[2]` → one of three call sites (`pickin.c:79`, `126-127`, or `167-168`, selected by the `option` argument) → `extra.s2`, the invariant mass$^2$ of the $(p_4,p_5)$ or $(p_3,p_4)$ subsystem depending on ordering choice.
- `rannums[3]` → `pickin.c:177`, `yy4 = cos(pi*rannums[3])` → feeds `levi.delta` at `pickin.c:184` (`levi.delta = delb-yy4*st*sqrt(dd)/(2*ap)`) → which fixes `extra.s1` at `pickin.c:185` (`extra.s1 = extra.t2+pickzz.w1+(2*dotp.p12*r3-4*levi.delta)/st`).

So the 4 inputs become the 4 invariants $t_1, t_2, s_2, s_1$ — matching
Appendix A of the Vermaseren paper. `pickin.c`'s own header comment
(`pickin.c:16-17`) names these directly:

```c
The variables used are in principle t1=(p1-p3)^2, t2=(p2-p5)^2,
s2=(p4+p5)^2 and s1=(p3+p4)^2.
```

and the same comment (`pickin.c:17-20`) is what motivates [Topic 3](topic3-subset-sum-proof.md) later:
it explicitly says the Jacobian is $1/\sqrt{-\Delta_4(p_1,p_2,p_3,p_4)}$
and that numerical stability improves by using $\Delta_4$ itself as an
integration variable instead of $s_1$ directly — this is the Jacobian-peak
"type of peak" flagged as an open thread under [Topic 1](topic1-peaks-and-O2n.md). **This is now
traced to its primary source, not just the lecture's citation list:**
Byers & Yang (1964, *Rev. Mod. Phys.* 36, 595), Theorem 13 (eq. 79),
derive the general-$n$ phase-space volume element in closed form as a
product of Gram-determinant-like invariants $\Delta_l$ raised to
$(n-5)/2$-type powers, exactly the object named $\Delta_4$ here — this
is the actual origin of "use $\Delta_4$ as an integration variable"
(`pickin.c`'s comment), matching the lecture's own history section
crediting Byers–Yang specifically for "solving phase space this way."

**Claim: the 5th ($3n-4=5$) variable — overall azimuth about the beam
axis — is computed, not sampled.**

`orient.c:68-69`:

```c
rr  = sqrt(-levi.gram/s)/(variab.p*variac.pp4);
variab.sp3 = rr/variac.pp3;
```

`sp3` (sin φ₃, per the file's own docstring at `orient.c:15`, "sp3 etc:
sim(phi3)") is computed algebraically from `levi.gram`, which in turn was
already fixed at `pickin.c:188` (`levi.gram = (1-yy4)*(1+yy4)*dd/ap`) purely
in terms of the 4 already-sampled invariants — no new random number enters
between `pickin` returning and `sp3` being set. The **only** remaining
freedom at this stage is a sign choice, resolved at `orient.c:74-79`:

```c
a1  = variac.pp3*variab.cp3-variac.pp5*variab.cp5;
if ( fabs(variac.pp4+variac.pp3*variab.cp3+variab.cp5*variac.pp5)
        >= fabs(fabs(a1)-variac.pp4) ) {
    if ( a1 <  0 ) variab.cp5 = -variab.cp5;
    if ( a1 >= 0 ) variab.cp3 = -variab.cp3;
}
```

This is a **deterministic consistency check** (which branch to take is
decided by comparing two already-computed magnitudes), not a random draw —
confirming there is no hidden 5th `rannums[]` consumption disguised here.

**Conclusion, now grounded in the quoted lines above rather than asserted:**
the code fixes the trivial 5th variable (event orientation about the beam
axis) to a convention — implicitly, φ₃ and φ₅ are pinned relative to each
other and to a plane, per the file's own remark that particle 4 determines
the reference plane (matches the lecture's aside about $e_{p_1p_2p_3p_4}$
and choosing p4 to fix the xz-plane) — rather than sampling and later
averaging over it, because doing so would only multiply the total by the
flat factor $2\pi$ for unpolarized beams. This is the concrete mechanism
behind "at most" in "at most $3n-4$": the code uses 4 of the 5 formally
available slots, having recognized the 5th carries no information worth
paying sampling cost for.

### Why "at most" (not "exactly") $3n-4$

The phrase "at most $3n-4$" (not "exactly") now reads as deliberate too:
$3n-4$ is the *ceiling* set by naive dimensional counting of the
phase-space integral, but a specific diagram/parametrization, once you
exploit its own symmetries (as `orient.c` does with the azimuthal angle
here), may need strictly fewer *sampled* variables than that ceiling — the
"missing" ones are accounted for by direct construction/convention rather
than Monte Carlo integration.

### Open threads / not yet resolved

- **Resolved** (previously open): whether the $n=3$ "one azimuthal
  direction is always free" finding generalizes to any $n$, not just our
  case. Byckling & Kajantie, *Particle Kinematics* (1973), Table III.1
  (p.49) states the general-$n$ result directly: "Final State Variables"
  ($3n-4$, includes one trivial azimuth $\phi$) reduces to "Essential
  Final State Variables" ($3n-5$, $\phi$ dropped) for *any* $2\to n$
  process — not a special feature of $n=3$. See the table added above.
  What remains genuinely open, and is *not* resolved by this table: for
  a $2\to4$ process built from a 2-body sub-decay (`eemm.c`/`gamgam.c`/
  `epmm.c`, where `rannums[4..6]` appear), the *additional* azimuthal
  random draw `phicm6` in `gamgam.c` is explicitly flagged in this
  document (Topic 2 caveats above) as **not** trivial in the same way —
  the sub-decay plane's orientation relative to the produced pair is
  physical once polarization/spin correlations matter. Byckling–Kajantie's
  table counts exactly *one* trivial azimuth per $2\to n$ process (the
  overall azimuth about the beam axis); it does not by itself say
  anything about additional azimuths introduced by sub-decay structure,
  which is a separate physical question this document has not resolved
  and would be worth its own topic if the $2\to4$ exercise is built.
