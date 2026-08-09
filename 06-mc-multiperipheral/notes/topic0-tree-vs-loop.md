# Topic 0: what "tree diagram" and "loop" actually mean

[Topic 3](topic3-subset-sum-proof.md) proves a fact about internal-line momenta in tree
diagrams, but uses "tree diagram" as if it were self-evident. This topic
was inserted *after* [Topic 3](topic3-subset-sum-proof.md) was written, once a question exposed that
"tree = no loops" needs its own careful definition, and that "no loops"
itself is not simply "no closed shape when you draw it," but a precise
statement about how many internal momenta are left undetermined. This
section works from first principles and corrects two mistakes made along
the way (both documented below, not hidden), because seeing what went
wrong is as instructive as the final definitions.

### Setting up the general framework

Every internal line of a Feynman diagram carries a momentum. Every vertex
enforces 4-momentum conservation: with a fixed reference direction assigned
to each line,

$$\sum_{\text{lines }\ell\text{ at vertex }v} \sigma_{\ell,v}\,q_\ell = 0,
\qquad \sigma_{\ell,v}=\begin{cases}+1 & \ell\text{ into }v\\-1&\ell\text{ out of }v\end{cases}$$

A diagram with $V$ vertices gives you $V$ such vector equations. But they
are never all independent: summing **all** of them together, every
internal line appears exactly twice (once from each of its two endpoints)
with opposite sign, so it cancels completely — what remains is exactly
overall momentum conservation of the external legs. That means summing
*all* $V$ vertex equations reproduces a fact you already know is true
(total incoming momentum = total outgoing momentum) rather than giving new
information. So there are only

$$V - 1$$

**independent** vertex constraints, not $V$. (This "one equation is always
redundant" fact is the same telescoping mechanism used in [Topic 3](topic3-subset-sum-proof.md)'s proof,
just applied to *all* vertices at once instead of one side of a single
cut.)

### A terminology warning: "external" is a graph role, not "physically fixed"

Before going further, one distinction has to be nailed down, because this
document uses "external momentum" in two senses that must not be
conflated, and an earlier draft of this document did conflate them
(flagged and fixed here): **"external" as used in this constraint-counting
argument, and in [Topic 3](topic3-subset-sum-proof.md)'s proof, means "external to the diagram graph"**
— a leg that is not itself an internal propagator solved for by the
vertex equations. It says nothing about whether that leg's momentum is
held fixed across an entire calculation.

Concretely, for our own diagram: $p_1,p_2$ (the beam momenta) *are*
physically fixed once and for all by the experiment (the beam energy).
But $p_3,p_4,p_5$ (the final-state momenta) are "external" only in the
graph-topological sense used here — they are not internal lines, so
[Topic 3](topic3-subset-sum-proof.md)'s proof treats them as the given data from which $q_1,q_2$ are
*derived*, for any *one* diagram evaluation. Across a full cross-section
calculation, however, $p_3,p_4,p_5$ are exactly what the phase-space
integral ranges over — they are not fixed at all once you zoom out from
"one diagram evaluation" to "the integral over all of phase space." Both
statements are true, about the same symbols, because they answer
different questions: "is $p_3$ solved-for by the vertex equations of
*this* diagram?" (no — it's external, in that sense) vs. "does $p_3$ take
one fixed numerical value throughout the calculation?" (no — it's an
integration variable). This distinction is used without further comment
from here on; see the "How can a tree diagram peak at all" subsection of
[Topic 1](topic1-peaks-and-O2n.md) for why it matters.

### The definition (constraint-counting, not graph-shape)

Let $E$ be the number of **internal** lines (propagators) in the diagram —
these are the unknowns, solved for by the vertex-conservation equations in
terms of whatever values the diagram's external legs happen to carry (see
the terminology warning just above — "given" here means "given once you
fix a point in phase space," not "fixed for the whole calculation").
Compare unknowns to independent constraints:

$$L \;=\; E - (V-1) \;=\; E - V + 1$$

- **$L=0$: every internal momentum is uniquely fixed** by the external
  momenta via the vertex equations alone. Nothing is left to choose or
  integrate over. **This is what "tree diagram" means.**
- **$L>0$: after using every independent constraint, $L$ internal momenta
  remain genuinely free** — not fixed by anything external. Quantum
  mechanically, *every* value such a free momentum could take contributes
  to the amplitude, so the Feynman rules require summing (integrating)
  over all of them: $\prod_{i=1}^{L}\int \frac{d^4k_i}{(2\pi)^4}$. **This
  is the actual definition of a loop, and $L$ is the loop order** — not
  "does the picture look like it has a closed shape," but "how many
  internal momenta does the vertex-constraint system leave undetermined."

The fact that $L$ also equals the number of independent cycles in the
diagram-as-a-graph is a **theorem**, not the definition — provable from
the constraint-counting picture above, and checked concretely below rather
than assumed.

### Correction #1 (made mid-discussion): a diagram is not "non-tree" just because $V>2$

An early wrong claim in this discussion was that "any diagram with more
than 2 vertices is automatically non-tree." This is false, and falsified
directly by explicit graphs (all checked with `networkx.is_tree`):

| diagram | $V$ | $E$ | $L=E-(V-1)$ | tree? |
|---|---|---|---|---|
| our multiperipheral diagram ($A,B,C$) | 3 | 2 | 0 | yes |
| 4-vertex **chain** $A$–$C$–$D$–$B$ | 4 | 3 | 0 | **yes** |
| 4-vertex **star** ($C$ hub, $A,B,D$ leaves) | 4 | 3 | 0 | **yes** |
| 4-vertex **box** (rectangle, see Topic "loop example" below) | 4 | 4 | 1 | no |
| 10-vertex random tree | 10 | 9 | 0 | yes |

Two of the rows have **the same vertex count ($V=4$)** yet different
answers (tree vs. not) — so $V$ alone cannot possibly be the deciding
factor; a 10-vertex diagram can be a perfectly good tree. What actually
decides it is $E$ **relative to** $V-1$: a tree always has exactly one
fewer internal line than it has vertices ($E=V-1$); anything with one more
internal line than that ($E=V$) has exactly one loop, and so on. This
matters because "just count vertices" was floated as a shortcut — it
isn't one; it gives the wrong answer on two of the five rows above.

### Correction #2 (made mid-discussion): a naive symbolic check can give a false negative

While verifying the tree case symbolically (below), an initial check of
"is vertex $C$'s equation implied by $A$'s and $B$'s?" returned **False**
in `sympy`. This was not a flaw in the argument — it was a bug in the
check: $p_4$ was left as an unconstrained free symbol instead of also
being expressed via overall conservation
($p_1+p_2=p_3+p_4+p_5$) before comparing. Once that substitution is made,
the check correctly returns `0` (true). Recorded here because it's a
useful lesson for any future symbolic verification in this repo: **always
impose overall external-momentum conservation explicitly before checking
whether a vertex equation is "automatically satisfied" — otherwise the
external momenta look independent when they are not**, and the check will
spuriously fail. The same bug appeared a second time (below, in the box
case, where it produced `EmptySet` instead of a clean "no solution" or
"one free parameter" answer) before being caught and fixed identically.

### Worked example 1: our diagram is a tree ($L=0$)

<img src="../figures/multiperipheral.svg" alt="multiperipheral tree with q1, q2 momentum-cut demonstration" width="430">

(source: `../figures/multiperipheral.tex`, `tikz-feynman` — same figure used
again in [Topic 3](topic3-subset-sum-proof.md)'s diagram section below, since it's the same diagram.)

Vertex equations (incoming-positive convention):

- Vertex $A$: $p_1 = p_3+q_1 \Rightarrow q_1=p_1-p_3$
- Vertex $B$: $p_2 = p_5+q_2 \Rightarrow q_2=p_2-p_5$
- Vertex $C$: $q_1+q_2=p_4$

$V=3$, so independent constraints $=V-1=2$. Solving using only vertices
$A,B$ (2 equations, 2 unknowns $q_1,q_2$) gives a **unique** solution,
matching `pickin.c:16`'s own `t1=(p1-p3)^2, t2=(p2-p5)^2`. Vertex $C$'s
equation, checked with overall conservation $p_1+p_2=p_3+p_4+p_5$ imposed:

```python
import sympy as sp
p1,p2,p3,p4,p5,q1,q2 = sp.symbols('p1 p2 p3 p4 p5 q1 q2')
sol = {q1: p1-p3, q2: p2-p5}                          # from vertices A,B
p4_expr = sp.solve(sp.Eq(p1+p2, p3+p4+p5), p4)[0]       # overall conservation
sp.simplify((q1+q2).subs(sol) - p4_expr)                 # -> 0: vertex C is automatic
```

$E=2$, $L=E-(V-1)=2-2=0$. **Tree**, confirming [Topic 3](topic3-subset-sum-proof.md)'s proof applies.

### Worked example 2: the box diagram is one loop ($L=1$)

Two fermion lines ($a_1\to a_2$ and $b_1\to b_2$) exchanging two photons,
forming a closed rectangle — a standard one-loop QED box (two-photon
exchange):

<img src="../figures/box_loop.svg" alt="box diagram, one loop" width="600">

(source: `../figures/box_loop.tex`, a repo-native `tikz-feynman` diagram —
originally sketched from a photo the user shared outside the repo, rebuilt
here rather than referencing that external file, using the standard
particle-physics drawing convention: straight fermion lines with arrows,
sine-wave photon propagators, compiled with `lualatex` and converted to
SVG.)

Two fermion lines run through the diagram: $a_1\to TL\to TR\to a_2$ (top,
labeled $k_1$ on the internal segment) and $b_1\to BL\to BR\to b_2$
(bottom, $k_2$), exchanging two photons $k_3$ (left) and $k_4$ (right).
The prose and algebra below use position-based names ($TL,TR,BL,BR$ for
the vertices; $k_{top}=k_1$, $k_{bottom}=k_2$ for the two fermion
segments, $k_{left}=k_3$, $k_{right}=k_4$ for the two photons) since those
map directly onto the proof; the figure labels are the same four internal
lines under their $k_1$–$k_4$ names. The loop
**is** the full rectangle, walked $TL\to TR\to BR\to BL\to TL$: cutting
any one of the four internal lines leaves the other three still
connecting all four vertices, so no single cut disconnects the diagram —
the graph-theoretic signature of $L\geq1$.

$V=4$, $E=4$ (one more than the $V-1=3$ a tree would have). Solving the
vertex-conservation system explicitly (imposing overall external
conservation $a_1+b_1=a_2+b_2$ first, per Correction #2 above):

```python
import sympy as sp
kt,kr,kb,kl = sp.symbols('k_top k_right k_bottom k_left')
a1,a2,b1,b2 = sp.symbols('a1 a2 b1 b2')
eqTL = sp.Eq(a1, kt+kl); eqTR = sp.Eq(kt+kr, a2)
eqBL = sp.Eq(b1+kl, kb); eqBR = sp.Eq(kb, b2+kr)
a2_expr = sp.solve(sp.Eq(a1+b1, a2+b2), a2)[0]
sp.linsolve([eqTL, eqTR.subs(a2,a2_expr), eqBL, eqBR], [kt,kr,kb,kl])
# -> {(a1 - k_left,  b1 - b2 + k_left,  b1 + k_left,  k_left)}
```

Every internal momentum ($k_{top},k_{right},k_{bottom}$) comes out
expressed *in terms of* `k_left` — but `k_left` itself is never pinned to
a value by any combination of the four vertex equations. That is the loop
momentum: exactly the free parameter the definition predicts, matching
$L=E-(V-1)=4-3=1$. Physically you'd write the amplitude with
$\int d^4k_{\text{left}}/(2\pi)^4$ left in it — a genuine one-loop QED box
diagram (two-photon exchange correction).

### Worked example 3: a 4-vertex diagram that is *still* a tree ($L=0$)

Direct answer to "why bother with $E-V+1$, why not just count vertices" —
here is a diagram with the *same* vertex count as the box above ($V=4$)
that is nonetheless a tree, because its 3 internal lines form a **chain**,
not a cycle:

<img src="../figures/chain4.svg" alt="4-vertex chain tree diagram" width="750">

(source: `../figures/chain4.tex`, `tikz-feynman`.) Process: $\gamma + e^- \to e^- + \gamma +
\gamma + \gamma$ — an incoming photon Compton-scatters off an electron,
which radiates two extra bremsstrahlung photons, drawn as one continuous
fermion line running through 4 vertices $A,C,D,B$, each vertex emitting
one photon. **Every vertex was checked to be a genuine QED vertex** (2
fermion legs + 1 photon leg each — an earlier draft of this diagram
mistakenly gave vertices $A$ and $B$ three fermion legs and zero photons,
which is not a valid QED vertex; fixed by making one of each end-vertex's
external legs a photon instead of a fermion):

| vertex | legs | valid? |
|---|---|---|
| $A$ | $e^-_{\text{in}}$, $\gamma_{\text{in}}$, $q_1$ (fermion, out) | 2 fermion + 1 photon — yes |
| $C$ | $q_1$ (fermion, in), $q_2$ (fermion, out), $\gamma_1$ (out) | 2 fermion + 1 photon — yes |
| $D$ | $q_2$ (fermion, in), $q_3$ (fermion, out), $\gamma_2$ (out) | 2 fermion + 1 photon — yes |
| $B$ | $q_3$ (fermion, in), $e^-_{\text{out}}$, $\gamma_3$ (out) | 2 fermion + 1 photon — yes |

Solving vertices $A,C,D$ (3 equations, 3 unknowns $q_1,q_2,q_3$) gives a
unique solution; vertex $B$'s equation is then automatically satisfied
once overall conservation is imposed (checked, returns `0`, same pattern
as worked example 1). $V=4$, $E=3=V-1$, so $L=E-(V-1)=3-3=0$: **tree**,
despite having the same vertex count as the one-loop box in example 2.

**This is the direct resolution of the "why not just count vertices"
objection:** examples 2 and 3 both have $V=4$; one is a tree, one is a
one-loop diagram. Vertex count cannot distinguish them — only the
relationship between $E$ and $V-1$ can, because that relationship is
precisely asking "does the vertex-conservation system have a unique
solution, or is something left over to integrate?", which is the actual
physical question at stake.

### Where the graph-theoretic "cycle" picture fits in

The geometric fact "cutting an edge from a tree always disconnects it,
cutting an edge from the box's rectangle does not" (visible directly by
inspection, and used informally earlier in this discussion) is not a
separate, alternative definition — it is a **theorem** that follows from
the constraint-counting definition above: a free/undetermined internal
momentum is precisely one you can route "around a closed path" without
violating any vertex equation, and closed paths in a connected graph exist
exactly when $E>V-1$. Both pictures were checked to agree in every example
above ($E-(V-1)$ arithmetic vs. `networkx.is_tree`/`cycle_basis`), which is
why either is a legitimate way to *check* tree-vs-loop in practice — but
the constraint-counting definition is the one that explains *why* it
matters physically (integration variables), and is the one to reach for
first when a diagram's shape is ambiguous or hard to draw.
