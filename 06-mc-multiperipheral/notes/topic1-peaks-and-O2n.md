# Topic 1: "peaks", "types of peaks", and why propagator count grows like $O(2^n)$

### The source sentence

> "...when there are many Feynman diagrams in a reaction, there can be many
> different types of peaks and the number of potential propagators that can
> cause peaking grows like $O(2^n)$, while the number of non-trivial
> integration variables is at most $3n-4$."
> — uam19 `part5.pdf`, p.1

### Parsing it precisely

The sentence makes **two separate claims**, and it's important not to
collapse them into one (an error corrected during this discussion):

1. **Qualitative claim:** "there can be many different **types** of peaks."
   This is about *kinds* of peaks, not a count.
2. **Quantitative claim:** "the number of potential **propagators** that can
   cause peaking grows like $O(2^n)$." The grammatical subject being counted
   is *propagators*; "that can cause peaking" is a restrictive clause
   narrowing which propagators are in view (only ones capable of going
   small/near-singular in the accessible kinematic region) — it is not a
   second, independently-counted quantity.

The two claims are related but distinct: claim 2 (propagator count) explains
*one* source of the "many different types" in claim 1, but not the only one.

### "Peak" vs. "pole"/"singularity" — deliberate, not sloppy

Vermaseren uses "peak" as the general/operative term and "pole" as a
narrower one — he writes "peaks/poles" later in the same paragraph, showing
he has both words available and chooses deliberately:

- A **pole** is a property of the integrand *function* — a point where a
  propagator denominator $q^2 - m^2$ hits exactly zero, making the
  integrand formally infinite.
- A **peak** is the *numerically relevant* property for Monte Carlo: a
  region where the integrand is large and sharply varying. Every pole
  produces a peak (if the accessible phase-space region gets close to it),
  but not every peak comes from a true pole — e.g. a Jacobian factor like
  $1/\sqrt{-\Delta_4}$ (see `orient.c`'s doc comment on the Gram
  determinant) peaks without any propagator going on-shell, and the
  numerator/denominator near-cancellation in `pi0.c` (formula 3.2, the
  $1/(t_1t_2)^2$ term) produces a numerically dangerous region without a
  literal pole there either.

So the "different types of peaks" (claim 1) plausibly include at least:

- **Propagator peaks** — $1/(q_i^2-m_i^2)$ small for some internal line $i$.
  This is what claim 2 counts, and what `mapt.c` (the `dt/t` mapping) and
  `t1`/`t2` in `pickin.c` exist to tame.
- **Phase-space-measure/Jacobian peaks** — e.g. the $1/\sqrt{-\Delta_4}$
  factor `orient.c` calls out, unrelated to any single propagator.
- **Cancellation peaks** — large near-equal terms nearly canceling (the
  $O(s^2 M_\pi^4)$ numerator terms in `pi0.c` that the whole Part 5 lecture
  is fundamentally about). Not a pole anywhere; numerically peak-like
  because relative precision collapses there.

### How can a tree diagram peak at all, if every momentum is uniquely determined?

A sharp objection worth resolving before going further: [Topic 3](topic3-subset-sum-proof.md) proved
that in a tree diagram, every internal-line momentum is a *fixed*,
uniquely-determined function of the external momenta (a signed
subset-sum). If nothing is left free or undetermined the way a loop
momentum is, where could a "peak" or a "non-trivial integral" possibly
come from?

**Resolution: "uniquely determined" and "integrated over" are not in
tension, once "external momentum" is read in the sense fixed in [Topic 0](topic0-tree-vs-loop.md)'s
terminology warning** ("external to the diagram graph," not "physically
fixed for the whole calculation"). [Topic 3](topic3-subset-sum-proof.md)'s proof fixes *a single point*
$p_1,\dots,p_5$ in phase space and derives $q_1,q_2$ from *that one*
choice — true, and not in question; that is what "external, in the
graph-topological sense" buys you at that one point. But a cross section
is not the integrand evaluated at one kinematic point; it is
$$\sigma \propto \int d\Phi_n(p_3,\dots,p_{n+2}) \; |\mathcal{M}|^2,$$
an integral over *every* kinematically allowed choice of the final-state
momenta. $t_1,t_2,s_2,\cos\theta_4$ ([Topic 2](topic2-3n-minus-4.md)'s $3n-4=5$, minus the
trivial azimuth) are not extra freedoms *inside* one diagram's
kinematics — they parametrize *which point of the final-state phase
space you are currently at*. For each such point, $q_1=p_1-p_3$ is
indeed uniquely fixed, exactly as [Topic 3](topic3-subset-sum-proof.md) proved; the peak is not a peak
"in $q_1$" holding externals fixed, it is a peak **in the integrand as a
function of where in phase space you are evaluating it**, as $t_1$ itself
(one of the integration variables) sweeps close to its own on-shell value.

**This is directly visible in the code, not just an abstract argument.**
`eepi.c:14-20` — the function VEGAS calls once per Monte Carlo
sample — takes only two things as genuinely fixed input: the beam energy
`eemminput.sq` and the particle masses. It takes four **freshly-drawn
random numbers** `rannums[0..3]` (uniform on $[0,1]$) as its actual
per-call input, and *everything else, including $t_1,t_2$ and hence the
final-state momenta $p_3,p_4,p_5$, is derived inside `pickin`/`orient`
from those random numbers, fresh every call*:

```c
extra.t1 = mapt(t1min,t1max,&dt1,rannums[0],1);      // pickin.c:93
```

`t1min`,`t1max` (`pickin.c:90-92`) depend only on the fixed beam energy
and masses — **not** on `rannums`. So as `rannums[0]` sweeps $[0,1]$
across many VEGAS calls, `t1` sweeps continuously across its *entire*
physically allowed range $[t1min,t1max]$. For any *one* call, $t_1$ (and
hence $q_1=p_1-p_3$) is a single fixed number — but across the full
Monte Carlo run, $t_1$ takes on every value in that range, including
values arbitrarily close to its physical edge.

**Where the actual singularity sits, and why it's reachable.** `mapt.c`'s
own doc comment (`mapt.c:15`) states the mapping's purpose directly:

```c
/*
    Assumes dt/t.   t < 0
*/
```

$t_1$ is a spacelike momentum transfer, always $t_1<0$, and can approach
$t_1\to0^-$ within the physically allowed range — exactly where the
propagator denominator in `pi0.c`'s matrix element,
`part1 = -64*levi.gram/(tt*tt)` with `tt = extra.t1*extra.t2`
(`pi0.c:9`), blows up as $1/t_1^2$. `pickin.c:93` calls `mapt(...,
rannums[0], 1)` with `type=1`, i.e. the $dt/t$ log-mapping is what
actually runs — built for exactly this reason: a *flat* (uniform,
`type=0`) sampling of $t_1$ would place very few Monte Carlo samples near
$t_1\to0^-$ (a small sub-interval, linearly), even though that is
precisely where the integrand is largest and contributes most of the
cross section — the classic Monte Carlo variance problem. The $dt/t$ map
instead concentrates sample density logarithmically near $t_1\to0^-$, so
enough samples land where the integrand actually matters.

**So, precisely:** a tree diagram peaking is not a contradiction of
"every momentum is uniquely determined" — it is a restatement of it, one
level up. Because $q_1$ is a *fixed, continuous function of $t_1$*
($q_1^2=t_1$, trivially, by definition), and $t_1$ ranges continuously
over an interval that gets arbitrarily close to $q_1$'s own on-shell
point as part of the physical integration domain, the *integrand* — not
any individual momentum — develops a sharp peak as a function of where
in that domain you evaluate it. Determinism of $q_1$ given $t_1$ is
exactly what *transmits* the propagator's would-be pole at $t_1=0$ into
a peak in the phase-space integral over $t_1$; if $q_1$ were *not* a
fixed function of the integration variables, there would be no mechanism
for a peak in $t_1$-space to exist at all.

### Deriving the $O(2^n)$ propagator count

**Correction — scope of this derivation.** Vermaseren's sentence says
"the number of potential propagators that can cause peaking" with **no
restriction to tree diagrams anywhere in the source text** — "Feynman
diagrams in a reaction" is unrestricted. An earlier version of this
derivation silently substituted "in any *tree* diagram" for the source's
unrestricted "Feynman diagrams" and never flagged that as a narrowing.
That was a mistake: it is not Vermaseren's claim, it is a restricted
claim of ours, and it needed to be stated and justified as such, not
presented as if it were the same thing.

**Why the derivation below is nonetheless restricted to trees, stated
honestly:** the subset-sum lemma this count is built on ([Topic 3](topic3-subset-sum-proof.md)) is
proved *only* for tree diagrams — it relies on removing one internal line
splitting the diagram into exactly two pieces (the unique-path property
of trees), which fails once $L\geq1$: a loop diagram has $L$ internal
momenta left undetermined by vertex conservation ([Topic 0](topic0-tree-vs-loop.md)'s
$L=E-(V-1)$), so an internal line's momentum in a loop diagram is
generally an external-leg subset-sum **plus an integer combination of
the free loop momenta** — not a fixed subset-sum at all, and its "peak"
occurs at some point along a continuous loop-integration contour, not at
one fixed kinematic value the way a tree propagator's does. Counting
those is a genuinely different (harder) combinatorial problem — related
to counting distinct loop-diagram topologies at each loop order, which
itself grows factorially/combinatorially with the number of external legs
and internal loops, not via the simple subset-counting argument below.

**Update: we did attempt this count** (below), since it is worth
resolving rather than leaving open — the tree-only restriction turned
out to understate what generalizes, once done carefully.

### Non-tree generalization: loop-momentum routing and Landau singularities

**Part A — the routing lemma generalizes cleanly.** For a diagram with
$L$ independent loop momenta $k_1,\dots,k_L$ (fixed by picking a spanning
tree of the diagram graph; the $L$ edges *not* in the spanning tree, the
"chords," each carry one free loop momentum), the same telescoping-sum
argument as [Topic 3](topic3-subset-sum-proof.md)'s proof — now applied along the spanning tree, with
each chord contributing its own loop momentum at the point it closes the
loop — gives: **every internal line's momentum is a signed subset-sum of
the external momenta, plus an integer-coefficient combination of the $L$
loop momenta.** This was verified directly (not just argued) on the box
diagram (`../figures/box_loop.tex`, $L=1$) by reusing the vertex-conservation
solve already in [Topic 0](topic0-tree-vs-loop.md)'s Worked example 2 and checking each of the 4
internal lines' loop-momentum coefficient is an integer:

```python
# ../scripts/loop_momentum_routing.py — reuses NOTES.md [Topic 0](topic0-tree-vs-loop.md)'s sympy solve
# {(a1 - k_left,  b1 - b2 + k_left,  b1 + k_left,  k_left)} for (k_top,k_right,k_bottom,k_left)
# k1 (top):    external part = a1,      loop coeff = -1
# k4 (right):  external part = b1 - b2, loop coeff =  1
# k2 (bottom): external part = b1,      loop coeff =  1
# k3 (left):   external part = 0,       loop coeff =  1   (k3 IS the loop momentum)
```

All four coefficients are integers ($\pm1$), confirming the claim on this
example. This part of the generalization is straightforward — it is
standard loop-momentum routing, not new physics.

**Part B — "how many can peak" does *not* generalize the same way, and
here is why, precisely.** A tree propagator's momentum is a *fixed*
function of the external momenta: for given external kinematics, it
either sits near on-shell or it doesn't — a discrete, checkable condition
per line, which is exactly what Topic 1's subset count enumerates. A loop
propagator's momentum depends on the loop-integration variable, which
ranges over all of (complexified) $\mathbb{R}^4$ per loop — so its
on-shell condition $q_e^2=m_e^2$ defines a whole *surface* in
loop-momentum space that is, generically, crossed by the integration
contour for **any** external kinematics. A single on-shell loop
propagator is therefore not by itself a singularity of the integral: the
$+i\epsilon$ (Feynman) prescription deforms the contour around it without
obstruction. This was checked explicitly, not asserted: solving
$q_{k_1}=0$ alone for the box's loop momentum $l$ gives $l=a_1$ for *any*
$a_1$ — no constraint on external kinematics at all.

A **genuine** singularity of the loop integral (a "Landau singularity")
requires the contour to be *pinched* — trapped between two or more
poles approaching from opposite sides — which needs **at least two**
on-shell lines closing off the one-complex-dimensional contour
simultaneously, and more generally **at least $L+1$ on-shell lines** to
pinch all $L$ loop-integration contours at $L$ loops (standard result;
see e.g. Eden–Landshoff–Olive–Polkinghorne, *The Analytic S-Matrix*).
This was verified directly for the box ($L=1$): solving $q_{k_1}=0$ and
$q_{k_3}=0$ *simultaneously* has no solution unless $a_1=0$ — a genuine,
nontrivial condition on external kinematics, unlike the single-line case.

**Full classification, computed and checked, not merely argued
(`../scripts/landau_singularity_box.py`):** for every one of the $2^4-1=15$
nonempty subsets of the box's 4 internal lines, solve for the loop
momentum from the on-shell conditions and check whether a nontrivial
condition on external kinematics survives:

```
S={k1}, {k2}, {k3}, {k4}   (|S|=1, all 4): generic, no kinematic
                                             condition -- NOT genuine singularities
S={k1,k2}: a1+b1=0        S={k1,k3}: a1=0          S={k1,k4}: a1+b1-b2=0
S={k2,k3}: b1=0           S={k2,k4}: b2=0          S={k3,k4}: b1-b2=0
                           (all 6 pairs: genuine codim-1 loci)
S of size 3 (all 4 triples): genuine codim-2 loci (intersections of two pair-conditions)
S={k1,k2,k3,k4}: genuine codim-3 locus (intersection of three pair-conditions)
```

11 genuine singularity loci out of 15 subsets — exactly matching the
closed-form prediction $\sum_{k=L+1}^{E}\binom{E}{k} = \sum_{k=2}^{4}\binom{4}{k}
= 6+4+1=11$ for $E=4,L=1$. The formula was checked against the explicit
enumeration by the script itself (`assert genuine == predicted`), not
just quoted.

**Sanity check — this correctly reduces to the tree case.** At $L=0$,
the rule "$|S|\geq L+1$" becomes "$|S|\geq1$": *every* single internal
line's on-shell condition already counts, with no floor excluded — which
is exactly Topic 1's tree picture, where a single propagator's on-shell
condition ($t_1=(p_1-p_3)^2=m^2$, say) is *directly* a constraint on
external kinematics, because there is no loop momentum to absorb it.
This consistency check passing (rather than being assumed) is what makes
the generalization trustworthy.

**The resulting general statement.** For a diagram with $E$ internal
lines and $L$ loops, the number of genuine (leading) Landau-singularity
loci is bounded by
$$\sum_{k=L+1}^{E}\binom{E}{k} \;\leq\; 2^E,$$
**the same $O(2^E)$ order of growth as the tree case** (which is the
$L=0,\,|S|\geq1$ special case, summing the *full* binomial range) — just
with the smallest-subset "floor" pushed up by $L$. So: **for a single
fixed diagram, going to higher loop order does not blow the exponential
past $O(2^E)$ in that diagram's own internal-line count** — it stays the
same order, only removing the (few, small) subsets that loop freedom
absorbs.

**What is genuinely new, and left open.** The above bounds singularities
*for one fixed diagram topology*. Vermaseren's $O(2^n)$, as derived for
trees in Topic 1, summed over subset-sums *across every diagram topology
contributing to the reaction* — and the number of *distinct diagram
topologies themselves*, at fixed external legs but increasing loop order,
grows combinatorially/super-exponentially in $L$ (a much harder counting
problem than anything solved here — related to the enumeration of
1-particle-irreducible graphs, not attempted in this document). So the
loop-momentum-routing and Landau-singularity results above answer "does
$O(2^n)$-type behaviour still hold for a single loop diagram" (yes, same
order, verified) but do **not** answer "what is the total peak count
summed over all diagrams at all loop orders contributing to a reaction"
— that remains open, and is plausibly a much larger quantity than
$O(2^n)$, not a smaller one. Given that (per the scope note above)
nothing in this course actually performs loop integration, this is
recorded as a resolved side-investigation, not something this exercise
needs to complete further.

**Why restricting to trees is a defensible reading of Vermaseren's intent
here, even though he does not say so explicitly:** every routine and
worked example in this lecture (`pickin`/`orient`/`eepi`, the
double-t-channel diagram, the "12 diagrams" aside for
$e^-e^+\to e^-e^+\mu^-\mu^+$ at "only QED") is tree-level phase-space
kinematics for Born-order cross sections — nothing in `part5.pdf` or the
surrounding `kinc2.tar.gz` archive performs loop integration. So in the
*practical context this lecture addresses*, "Feynman diagrams in a
reaction" most likely does mean the set of tree diagrams contributing at
that order — but this is our inference from the lecture's scope, not
something the source text states, and it should not be read as
Vermaseren having limited his own claim; he simply never needed to
address loop diagrams because none appear anywhere in this material.

**Setup.** For an $n$-particle final state produced from a $2\to n$ process,
there are $n_{ext} = n+2$ external legs (2 incoming + $n$ outgoing).

**Key fact (momentum conservation on a tree diagram — see the scope note
above for why this derivation is restricted to trees):** the momentum
flowing through *any* internal line, in *any tree diagram* for this
process, equals (up to overall sign) the sum of a subset $S$ of the
external momenta — this follows just from conservation at each vertex,
regardless of which specific tree topology you draw. A subset $S$ and its
complement $S^c$ (all other legs) describe the *same* physical line,
since $q$ and $-q$ give the same propagator denominator $q^2-m^2$.

**Counting distinct possible propagator momenta.** The number of subsets of
$n_{ext}$ legs is $2^{n_{ext}}$. Divide by 2 for the $S\leftrightarrow S^c$
symmetry, then exclude the trivial subsets that aren't internal lines at all
(the empty set / full set — that's just overall momentum conservation — and
single-leg subsets / their complements, which are external legs, not
internal propagators). Working this out exactly:

$$\text{distinct propagator momenta} = 2^{n_{ext}-1} - n_{ext} - 1 = 2^{n+1} - n - 3$$

This was verified by brute-force enumeration (see calc below) for
$n=1,\dots,7$ — exact match. Asymptotically this is $\sim 2^{n+1}$, i.e.
$O(2^n)$: Big-O absorbs the factor of 2 and the $-n-3$ correction, which is
the sense in which this *tree-restricted* count matches the lecture's
$O(2^n)$ — an order-of-growth statement, not an exact formula, and (per
the scope note above) a lower bound on Vermaseren's unrestricted claim,
not a proof of it: loop-diagram propagators are not counted here at all.

**A gap this closes, not previously checked:** the count above is a
count of abstract bipartitions of $n_{ext}$ legs — it was never verified
that every such bipartition is actually *achievable* as an internal line
of some real cubic (trivalent) tree topology, as opposed to including
bipartitions no valid diagram could ever realize (which would make the
count an overcount). Checked directly for $n=3$
(`../scripts/propagator_achievability_check.py`): enumerating every cubic
tree topology on the 5 external legs and collecting every bipartition
realized by *some* internal line of *some* topology gives exactly the
same 10 bipartitions the formula counts — confirmed by exact set
equality, not just matching cardinality. So the formula is not an
overcount; every counted propagator momentum genuinely occurs in at
least one tree diagram for the process.

| $n$ (final-state particles) | $n_{ext}=n+2$ | exact distinct propagator momenta | $2^n$ |
|---|---|---|---|
| 1 | 3 | 0 | 2 |
| 2 | 4 | 3 | 4 |
| 3 | 5 | 10 | 8 |
| 4 | 6 | 25 | 16 |
| 5 | 7 | 56 | 32 |
| 6 | 8 | 119 | 64 |
| 7 | 9 | 246 | 128 |

**What $n$ is:** the number of final-state particles — confirmed directly
by the very next clause in the same sentence, "the number of non-trivial
integration variables is at most $3n-4$", which is the standard
Byckling–Kajantie phase-space dimension count for an $n$-body final state
(sanity check: $n=2 \Rightarrow 3(2)-4=2$, the two independent angles of a
2-body decay in its rest frame — correct).

### Why this matters: the actual point of the sentence

The number of **potential** peak-causing propagators grows *exponentially*
in $n$ ($O(2^n)$), while your budget of **integration variables** to build
importance-sampling mappings along grows only *linearly* in $n$ ($3n-4$).
That mismatch is the structural reason multiperipheral/many-particle phase
space is hard: you cannot in general build one universal coordinate system
that flattens every potential peak using only $3n-4$ coordinates once the
number of potential peaks outstrips that budget. Hence the very next
sentence in the lecture: split phase space into pieces, each handling only
a *limited* number of peaks/poles, with a dedicated parametrization per
piece.

### Why splitting phase space "resolves" the $O(2^n)$ vs. $3n-4$ mismatch

**Correction to an earlier draft of this section:** an earlier version of
this note argued that $O(2^n)$ "has nothing to do with algorithmic time
complexity" and was purely a combinatorial/structural count. That was
wrong, and the lecture's own opening sentence says so directly: "this
does *not* mean that with the routines we study here we can calculate all
reactions trivially" — stated as the immediate consequence of the
$O(2^n)$-vs-$3n-4$ mismatch. This is explicitly a claim about
computational tractability, not just a taxonomy of peak types. The
correct relationship is: $O(2^n)$ *is* directly a time-complexity-driving
quantity, via the following chain, not a separate concern from it.

**Scope reminder, carried forward from "Deriving the $O(2^n)$ propagator
count" above, and easy to lose track of in what follows:** everything
below uses "$O(2^n)$" loosely for the *tree-restricted* count
$2^{n+1}-n-3$ we actually derived and verified. That count was explicitly
established as a **lower bound** on Vermaseren's unrestricted claim, not
a proof of it — loop-diagram contributions are not included. Nothing
below should be read as having tightened that to an equality; the
argument's *qualitative* structure (exponential hazard count vs. linear
coordinate budget) holds either way, but the specific number used for
"$O(2^n)$" in the concrete examples below (e.g. "10 potential propagators
for $n=3$") is the tree-only figure, not a verified count of all
diagrams' propagators.

**The chain from "potential peaks" to "computation time."** Each of the
$O(2^n)$ potential propagator momenta is a location where, if the
diagram/kinematic region makes that propagator go near on-shell, the
integrand becomes sharply peaked. Plain (unweighted/uniform) Monte Carlo
sampling converges as $O(1/\sqrt{N})$ *only if the sampling density
roughly tracks the integrand* — near an untamed peak it does not, so the
variance of the estimator blows up and the number of samples $N$ needed
to reach a fixed target accuracy grows correspondingly (in the worst
case, without any importance-sampling map at all, an integrable but sharp
peak can require many orders of magnitude more samples, or fail to
converge to acceptable precision in practical time at all). Importance
sampling (`mapt.c`'s `dt/t`, `mapw.c`'s `dw^2/w^2`, etc. — [Topic 2](topic2-3n-minus-4.md)'s code
citations) is precisely the fix, but each such map consumes one of your
$3n-4$ integration-variable "slots" to concentrate sampling density along
one specific singular direction. So: $O(2^n)$ counts the number of
distinct places variance could blow up, $3n-4$ counts how many of those
places you can afford to build a dedicated map for *in one pass*. The
mismatch between an exponentially-growing hazard count and a
linearly-growing map budget is exactly why "calculate all reactions
trivially" fails as the number of diagrams grows — this is a direct,
not incidental, algorithmic-complexity statement.

**What $O(2^n)$ counts, re-examined.** Go back to Topic 1's own
derivation: $2^{n+1}-n-3$ is a count of *distinct possible propagator
momenta across every diagram/topology that could contribute to the
reaction* — a property of the reaction ($n$ final-state particles), not
of any single diagram or any single Monte Carlo run. The worked example
made this concrete: for $n=3$ there are 10 potential propagator momenta
in total, but the *one* double-t-channel diagram this lecture builds
kinematics for uses only 2 of them ($t_1,t_2$). Critically, this does not
make $O(2^n)$ a purely combinatorial/non-computational fact: *every one*
of those 10 potential propagators is a place some diagram contributing to
the full reaction could peak, and a complete calculation of the reaction
(summing all contributing diagrams) must, in some piece of the
calculation, be numerically safe near each of them. $O(2^n)$ is a
statement about the size of the *menu*, and the size of that menu is
exactly what drives total computation cost once you account for every
diagram, not just the one worked example covers.

**Where the real constraint bites — an assumption now backed by the
primary source, not just illustrated by one example.** A single Monte
Carlo integration — one fixed set of $3n-4$ coordinates, one
importance-sampling map built from those coordinates — tames
singularities by dedicating one coordinate's map to one propagator
(exactly what `mapt.c`'s `dt/t` mapping does for $t_1$, or `mapw.c`'s
`dw^2/w^2` for an $s$-channel invariant — see the [Topic 2](topic2-3n-minus-4.md) code
citations above, and note `pickin.c` does exactly this: one coordinate
each for $t_1,t_2,s_2$). The argument below relies on this being
essentially one-coordinate-per-peak. An earlier version of this document
flagged that as an unproven assumption, illustrated only by the
$t_1,t_2$ example. Having now read Vermaseren's paper directly (p.349;
quoted in full in the "What splitting actually buys you" section below),
this is confirmed as the actual design of "event generation of the third
kind" — the method `pickin`/`orient` implements: the propagator
denominators are *directly* used as integration variables, one each, and
mapped via $\ln(-t_1),\ln(-t_2)$, handling both peaks together in one
unsplit calculation precisely *because* each gets its own coordinate.
The paper does not claim this is the only possible construction (a
cleverer joint map might exist in principle), but it is the method
actually used and documented, not merely an assumption read into the
code. Taking it as the operative design: you have $3n-4$ coordinates and
therefore can align with, at best, on the order of $3n-4$ independent
singular directions in one such construction. Once the number of
*potential* peaks in play exceeds that (which $O(2^n)$ guarantees
happens for any reaction with enough diagrams), no single coordinate
system can be tangent to all of them — some peaks are necessarily left
unmapped, and the integrand blows up in a direction your sampling density
does not follow, which is exactly the numerically catastrophic scenario
`pi0.c`'s $1/(t_1t_2)^2$ term represents (Topic 1's "cancellation peak"
case) if $t_1,t_2$ are not handled by dedicated maps.

**Correction, now grounded in the primary source rather than inferred:**
an earlier version of this section claimed `pickin`/`orient` itself *is*
an example of "splitting phase space into pieces." Having now read
Vermaseren's own paper (J.A.M. Vermaseren, "Two-photon processes at very
high energies," *Nucl. Phys.* B229 (1983) 347–371, p.349 — the very paper
Vermaseren's `part5.pdf` lecture cites and this whole exercise is based
on), that was wrong: the paper explicitly describes **two different
methods**, and `pickin`/`orient` is a worked example of the *other* one,
not the splitting method.

- **"Event generation of the second kind"** (p.349, first method):
  rewrite the phase-space integral so experimental cuts each become a
  single integration variable. Disadvantage, in the paper's own words:
  "the peaks of the matrix element now show up as a correlation between
  several variables so it is harder to integrate over them. It might be
  necessary therefore **to split the phase space up into various pieces
  such that each piece contains a single peak**. One can then
  concentrate separately on each piece in order to integrate it
  accurately." This is *the* splitting method Topic 1's source sentence
  refers to — and the paper is explicit that pieces are organized
  **one peak per piece**, confirming (from the primary source, not
  inference) the open question left in the previous version of this
  document about what a "piece" is organized around.
- **"Event generation of the third kind"** (p.349, second method — **the
  one the paper actually uses, and what `pickin`/`orient` implements**):
  rewrite the integral so the propagator denominators $t_1,t_2$
  themselves become the integration variables, then map them via
  $\ln(-t_1),\ln(-t_2)$ — exactly `mapt.c`'s `dt/t`. Combined with VEGAS's
  own adaptive multi-dimensional grid, this handles **both** $t_1,t_2$
  peaks *in one unsplit calculation*, without partitioning phase space at
  all. The paper states its advantage directly: "the great advantage of
  this method is its universality as one program can deal with any kind
  of experimental cuts" — the opposite tradeoff from splitting, which the
  paper says needs "a completely different reformulation of the
  phase-space integrals if the nature of the cuts is changed drastically"
  (a disadvantage explicitly *not* shared by the method `pickin`/`orient`
  uses).

So `pickin`/`orient` is not a "dedicated piece" of a split calculation —
it is a complete, unsplit generator for the one double-t-channel diagram,
using the *other* of the two methods the paper describes specifically to
avoid needing to split at all for this diagram's two peaks. This also
sharpens (rather than merely leaves open) the earlier flagged
"one-coordinate-per-peak" assumption: `pickin`/`orient` is direct
evidence that **within the $\ln(-t)$-style method**, one dedicated
coordinate per propagator peak, handled together in a single VEGAS
integration, is not just possible but is exactly the paper's preferred
approach for a small, fixed number of peaks (here, two) — splitting
(method one, single-peak-per-piece) is reserved for the harder case where
cuts *correlate* several variables' peaks together, not simply "too many
peaks for the coordinate budget" as a previous version of this document
suggested. Whether the true multi-diagram case (many diagrams, each
contributing its own $t$-type peaks) is better handled by generalizing
`pickin`/`orient`'s method (more coordinates, still unsplit) or by
splitting into single-peak pieces is a real design choice the paper does
not resolve in general — it resolves it only for this one diagram.

**Correction — a citation this document previously got wrong.** An
earlier version of this section attributed the splitting-into-pieces idea
itself to "Byers–Yang/Byckling–Kajantie," without having read either
source. Having now read Byers & Yang (1964, *Rev. Mod. Phys.* 36, 595 —
the paper cited in `part5.pdf`'s own history section for "solving phase
space this way"), that attribution was not supported: Byers–Yang is a
rigorous mathematical paper about (a) classifying which regions of
pairwise-invariant space $\{i,j\}=p_i\cdot p_j$ are physically realizable
by an actual set of $n$ four-momenta (via matrix rank/signature theorems;
the $n=3$ Dalitz-plot cubic surface, with its "horns" and conical points,
is the worked illustration), and (b) deriving the general-$n$
**phase-space volume element / Jacobian** in closed form (their Theorem
13, eq. 79) in terms of Gram-determinant-like invariants $\Delta_l$ —
**this** is what `pickin.c`'s own header comment is citing when it names
$\Delta_4$ and the Jacobian $1/\sqrt{-\Delta_4}$: Byers–Yang's formula is
the origin of using $\Delta_4$ as an integration variable for numerical
stability, not merely "credited in the lecture" as previously stated
without verification. Nowhere in Byers–Yang is there any discussion of
$O(2^n)$-style peak counting or splitting phase space into pieces — that
idea, per the actual text of Vermaseren's paper (p.349, quoted above), is
not attributed to a citation there either; it is presented as
Vermaseren's own methodological framing (his "event generation of the
second/third kind" terminology). The lecture's history bullet crediting
Byers–Yang for "solving phase space this way" refers specifically to the
invariant-variable/Jacobian formulation, not to the peak-counting
argument this document has been trying to trace to a citation.

**Where $\Delta_4$ and the Jacobian actually come from — Vermaseren's
Appendix A, formulas (A.1)-(A.3).** The paragraph above cites `pickin.c`'s
header comment and Byers-Yang's Theorem 13 for the origin of using
$\Delta_4$ as an integration variable, but had not quoted the specific
formulas `pickin.c` implements. uam19 `part5.pdf` names them directly:
*"Let us start with looking at the formulas A.1 and A.2. A.1 is the basic
formula..."* — these are in Appendix A of Vermaseren 1983 (p.360-361), not
in the lecture notes themselves. The two forms of the reformulated
$2\to3$ phase-space integral are:

$$\text{(A.1)}\qquad d\Phi = \frac{\pi}{16\,\lambda^{1/2}(s,m_1^2,m_2^2)}
\int\frac{ds_1\,ds_2\,dt_1\,dt_2}{\sqrt{-\Delta_4(p_1,p_2,p_3,p_4)}}$$

$$\text{(A.2)}\qquad d\Phi = \frac{\pi}{4\,\lambda^{1/2}(s,m_1^2,m_2^2)}
\int\frac{d\Delta\,ds_2\,dt_1\,dt_2}{(s_2+q-m_2^2)\sqrt{-\Delta_4(p_1,p_2,p_3,p_4)}}$$

(A.2) trades $s_1$ for the paper's "more exotic variable"
$\Delta\equiv p_1\!\cdot\!p_2\,q_1\!\cdot\!q_2-p_1\!\cdot\!q_2\,p_2\!\cdot\!q_1$
— this is exactly `levi.delta`, computed at `pickin.c:184`
(`levi.delta = delb-yy4*st*sqrt(dd)/(2*ap);`) and used at `pickin.c:185`
to obtain `extra.s1` from it, i.e. `pickin.c` actually integrates in the
(A.2) variable set ($\Delta,s_2,t_1,t_2$), then converts to $s_1$
algebraically rather than the other way around — consistent with the
paper's own framing of $\Delta$ as the numerically preferable choice.
Formula (A.3), immediately following, identifies $\Delta_4$ explicitly as
a Gram determinant and ties it to the Levi-Civita form used throughout
this document and [Topic 5](topic5-s-scaling-pi0.md):

$$\text{(A.3)}\qquad \Delta_4\equiv\varepsilon^{p_1p_2p_3p_4}\varepsilon_{p_1p_2p_3p_4}
=\varepsilon^{p_1q_1p_2q_2}\varepsilon_{p_1q_1p_2q_2}$$

— i.e. `pi0.c`'s `part1` numerator (`levi.gram`, formula (3.2)'s first
term, [Topic 5](topic5-s-scaling-pi0.md)'s subject) and `pickin.c`'s
phase-space Jacobian denominator (`orient.c`'s $1/\sqrt{-\Delta_4}$,
cited above as the canonical example of a Jacobian peak) are **the same
quantity**, $\Delta_4$, appearing in two different roles: as the
numerator of the matrix element's dominant term, and as the (inverse
square root of the) phase-space measure. This is not a coincidence
Vermaseren's paper treats as incidental — formula A.11
($D_1=\varepsilon^{p_1p_2q_2\mu}\varepsilon_{p_1p_2q_2\mu}$, p.361) notes
explicitly that $D_1$ (one of the boundary quantities `pickin.c` computes,
`levi.dd1`) "is a minor of the Gram determinant," the same
minor-of-a-determinant relationship formula (3.2)'s p.353 commentary uses
for `part2`/`part3` relative to `part1` (see [Topic
4](topic4-pseudoscalar-coupling.md)). Formulas (A.4)-(A.13), which follow
(A.3) in the paper, derive the $t_1,t_2,s_2$ boundary expressions
`pickin.c` implements almost line-for-line (e.g. `pickin.c:90-92`'s
`t1max`/`t1min` match (A.4)/(A.5)'s structure: a symmetric max expression
built from two $\lambda^{1/2}$ factors, and a numerically-stable min
expression built as a ratio to avoid the cancellation `pickin.c`'s own
header comment warns about). Working through the full correspondence
line-by-line for every boundary formula is not done here — flagged as an
open thread below — but the identification of (A.1)-(A.3) as the actual
source formulas, previously missing from this document despite citing
their consequences ($\Delta_4$, the Jacobian, `levi.delta`) repeatedly,
closes a real gap: earlier passages here referenced "the Gram determinant
$\Delta_4$" as if its definition were self-evident, without ever quoting
where in the primary source it is defined or how it is derived.

**Correction — "number of pieces," "number of peaks," and "number of
diagrams" are three different quantities, not one; an earlier draft of
this paragraph conflated them.** It claimed in one breath that the
number of pieces needed is "on the order of" $O(2^n)$ (the
propagator-subset count) and, two sentences later, that it's "tied to
$O(2^n)$... since a reaction with more diagrams needs more dedicated
configurations" — using "number of diagrams" and "number of potential
propagator momenta" as if they were the same growing quantity. They are
not, and the gap is large, not a rounding error: the number of distinct
cubic tree *topologies* on $n_{ext}$ labeled external legs is the
classical count $(2n_{ext}-5)!!$ (double factorial), which grows
**super-exponentially** — for $n=7$ ($n_{ext}=9$) it is $135135$,
already 550× larger than the propagator-subset count $246=2^{n+1}-n-3$
for the same $n$ (checked directly by enumeration, not just quoted).
Topic 1's own worked example already shows why these diverge: *one*
diagram (ours) realizes only $n-1=2$ of the 10 potential propagators for
$n=3$ — a diagram-count and a propagator-count are counting different
things (how many topologies exist vs. how many distinct lines could
appear across all of them), and nothing in this document establishes
they track each other.

**So, precisely, and now settled by the primary source rather than left
open:** splitting is a direct answer to an algorithmic-cost problem, not
a way of making that cost disappear. Vermaseren's paper (p.349, quoted
above) settles which of the two candidate quantities a "piece" is
organized around: **pieces are peak-organized, one peak per piece** —
"split the phase space up into various pieces such that each piece
contains a single peak," in the paper's own words — not diagram-organized.
So of the two counts distinguished above, (a) the number of potential
peaks (the tree-restricted $O(2^n)$ propagator-subset count, itself only
a lower bound on Vermaseren's unrestricted claim) is the one that
directly governs the number of pieces the *splitting* method (event
generation of the second kind) would need in the worst case; (b) the
number of diagrams is a related but separate driver of cost, relevant
because more diagrams generally means more potential peaks to organize
into pieces (each diagram contributes its own propagators to the menu),
and separately relevant to the *amplitude-vs-diagram-sum* tradeoff below
— but it is not itself "the number of pieces." The lecture's own aside on
p.3 (12 diagrams for $e^-e^+\to e^-e^+\mu^-\mu^+$, amplitude-squaring
preferred at that count) is evidence for a related but distinct
tradeoff — diagram count driving the cost of summing diagrams
*algebraically* before or after squaring — not direct evidence for how
many single-peak pieces a phase-space split would need for that
reaction; that would require knowing how many distinct peaks (in
Vermaseren's sense) those 12 diagrams' propagators produce, which this
document has not attempted to enumerate.

**A concrete instance in the archive of a *third* technique, not the
same as either of Vermaseren's two named methods — correcting an earlier
mislabeling.** An earlier version of this section called the `mgoto2`/
`gamgam` rotate-align-rotate-back trick (below) an instance of "dedicated
configuration" in the splitting sense. Having now read the paper's
precise two-method taxonomy above, that label doesn't fit either method:
it isn't splitting phase space into single-peak pieces (method one), and
it isn't rewriting the integral so a propagator denominator becomes the
integration variable directly (method two, what `pickin`/`orient` does).
It is a third, simpler technique — a coordinate transformation
(rotate/boost) applied *before* calling an otherwise-generic sampling
routine, so that routine's already-known blind spot lines up with the
physically relevant peak direction for this specific reaction. `mgoto2.c`
is a generic two-body decay routine, and its own doc comment
(`mgoto2.c:8-9`) states plainly:

```c
/*
    Routine for a two body decay. Note that there is no protection against
    special peaking behaviour near ct = +-1.
```

i.e. `mgoto2` is *not* itself built to handle any particular peak
direction — it is a general-purpose routine, left unmapped for the
$\cos\theta=\pm1$ danger. The lecture's own text (p.2) describes the fix
for the muon-pair case: the `gamgam` routine "rotate[s] the m4 system
first to align the photons along the z-axis, then do[es] the two body
decay and then rotate[s] and boost[s] back" — i.e. rotate the physically
relevant peak direction (incoming virtual photons) onto the coordinate
axis `mgoto2` already knows how to handle ($\cos\theta=\pm1$), sample,
then rotate/boost back. Neither `gamgam.c` nor `epmm.c` ship in
`kinc.tar.gz`/`kinc1.tar.gz`/`kinc2.tar.gz` (checked directly — absent
from all three archives), so this remains corroborating, not
directly-observed, evidence: `mgoto2.c`'s own docstring is the
first-hand code evidence; the rotate-align-rotate-back mechanism itself
is reported only in the lecture prose, not verified against a routine we
can read.

**One thing this is *not*:** `pickin.c`'s `option` parameter ([Topic 2](topic2-3n-minus-4.md)'s
code citations, `pickin.c:24-31`) is *not* an instance of any of the
three techniques above (splitting into single-peak pieces, the
$\ln(-t)$-propagator-as-variable method, or the rotate-align-rotate-back
trick), even though it superficially looks similar (a parameter that
changes which mapping function runs). It only changes the *order* of
integration ($s_2,t_1,t_2$ vs. $t_1,s_2,t_2$ vs. $t_1,t_2,s_2$) and which
map ($ds_2/s_2$ vs. a $\lambda$-function map) is used for $s_2$, always
for the *same* single double-t-channel diagram, within the one
$\ln(-t)$-method calculation. It tunes one detail of that one method's
own map; it does not select between different pieces, methods, or
diagrams. Conflating the two would be the same error flagged earlier in
this document (Correction #1, [Topic 0](topic0-tree-vs-loop.md)) — inferring more from a
superficial resemblance than the code actually supports.

### Worked example: our own diagram, $n=3$

Final state $p_3,p_4,p_5$ from $p_1,p_2 \to p_3,p_4,p_5$ (5 external legs).
Brute-force enumeration of all non-trivial leg subsets (mod
complementation) gives exactly **10** distinct possible propagator momenta
— matching $2^{n+1}-n-3 = 2^4-3-3=10$. Listed by which legs they group:

```
{p1,p2}  {p1,p3}  {p1,p4}  {p1,p5}
{p1,p2,p3}  {p1,p2,p4}  {p1,p2,p5}
{p1,p3,p4}  {p1,p3,p5}  {p1,p4,p5}
```

(each subset's complement is the same physical line, so only one
representative per pair is listed).

Of these 10 *potential* propagators, the **specific double-t-channel
diagram** this lecture builds kinematics for uses only **2**:

- $q_1 = p_1-p_3$, i.e. subset $\{p_1,p_3\}$ — this is `t1` in `pickin.c`
- $q_2 = p_2-p_5$, i.e. subset $\{p_2,p_5\}$ — this is `t2` in `pickin.c`

The other 8 subsets correspond to propagators that would appear in *other*,
undrawn diagrams for the same $e^-e^+\to p_3p_4p_5$ final state (e.g. an
s-channel-like diagram grouping $\{p_1,p_2\}$). This is the concrete
illustration of the gap between "propagators in one diagram" (small, grows
roughly linearly — $n-1=2$ here for a cubic tree) and "potential
propagators across all diagrams for the process" (the $O(2^n)$ ceiling the
lecture is bounding).

### Reproduction script

```python
from itertools import combinations

def distinct_propagator_momenta(n_ext):
    legs = list(range(n_ext))
    seen = set()
    for size in range(2, n_ext - 1):   # exclude 0,1,n_ext-1,n_ext (trivial)
        for combo in combinations(legs, size):
            s = frozenset(combo)
            comp = frozenset(legs) - s
            canon = min(s, comp, key=lambda x: sorted(x))
            seen.add(canon)
    return len(seen)

for n in range(1, 8):
    n_ext = n + 2
    print(n, n_ext, distinct_propagator_momenta(n_ext), 2**(n_ext-1) - n_ext - 1)
```

### Open threads / not yet resolved

- We have not yet catalogued the *other* "types" of peaks beyond the three
  listed above (propagator / Jacobian / cancellation) — is that an
  exhaustive list, or are there more categories the lecture has in mind
  later on?
- The $2^{n+1}-n-3$ count above is proved only for tree diagrams (per the
  scope note earlier in this Topic). **Update:** the "Non-tree
  generalization" subsection below now resolves the *per-diagram* part
  of this — a single loop diagram's genuine singularity count is still
  $O(2^E)$ in its own internal-line count $E$, verified via the Landau
  pinch classification on the box diagram. What remains genuinely open
  is the *sum over all diagram topologies at all loop orders* for a
  fixed reaction — the number of distinct topologies itself grows
  combinatorially with loop order, and we have not attempted that count.
  We also still have not confirmed that the tree-only reading is
  actually what Vermaseren has in mind rather than a genuinely
  unrestricted count over all orders — only argued it's reasonable given
  this lecture's entirely tree-level scope, not proved it's intended.
- An earlier draft of "Why splitting phase space resolves..." conflated
  "number of diagrams contributing to a reaction" with "number of
  potential propagator momenta" ($O(2^n)$ in the tree-restricted sense
  above) as if they were the same growing quantity driving the number of
  pieces a calculation needs. They are not: the number of distinct tree
  topologies on $n_{ext}$ legs grows like $(2n_{ext}-5)!!$
  (super-exponential — 549× the propagator-subset count already by
  $n=7$, verified in `../scripts/topology_count_vs_propagator_count.py`),
  not like $2^n$. This has been flagged and the conflation corrected
  (both are now presented as distinct, real concerns, not one
  quantity), but which one actually governs "how many dedicated
  phase-space pieces does a real calculation need" is not resolved here
  — it plausibly depends on the specific splitting/multichannel strategy
  used, which this document has not investigated.
