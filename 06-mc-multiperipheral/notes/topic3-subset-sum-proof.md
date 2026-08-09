# Topic 3: rigorous proof that an internal-line momentum is a signed subset-sum of external momenta

[Topic 1](topic1-peaks-and-O2n.md)'s $O(2^n)$ derivation rested on a claim stated but not proved: "the
momentum flowing through any internal line, in any tree diagram for this
process, equals (up to overall sign) the sum of a subset of the external
momenta." This section proves it rigorously and checks it against our own
diagram, both algebraically and graph-theoretically.

### Setup

A Feynman **tree diagram** is, as a graph, a tree: connected, no cycles.
Assign every line (internal or external) a momentum and a fixed reference
direction. At every vertex $v$, 4-momentum conservation reads

$$\sum_{\text{lines }\ell\text{ at }v} \sigma_{\ell,v}\, q_\ell = 0,
\qquad \sigma_{\ell,v} = \begin{cases}+1 & \ell\text{ points into }v\\ -1 & \ell\text{ points out of }v\end{cases}$$

### Claim

For any internal line $e$ of a tree diagram, removing $e$ splits the tree
into exactly two connected pieces $A,B$, which partition the external legs
into disjoint sets $S_A \sqcup S_B$. Then $q_e = \pm\sum_{i\in S_A} p_i$
for a sign fixed by the reference-direction convention.

### Proof

**Step 1 (pure graph theory — no physics yet).** A tree with $V$ vertices
has exactly $V{-}1$ edges and a *unique* path between any two vertices (two
distinct paths between the same pair would form a cycle, contradicting
"tree"). Removing one edge $e=(u,v)$ therefore leaves no path between $u$
and $v$, so the graph splits into exactly two connected components: $A$
(containing $u$) and $B$ (containing $v$). Every external leg attaches to
exactly one vertex, hence lies in exactly one of $A$, $B$ — this defines
the partition $S_A \sqcup S_B$.

**Step 2 (telescoping sum — the physics enters here).** Sum the
conservation equation over every vertex $v \in A$. Any internal line with
*both* endpoints inside $A$ contributes to exactly two of the summed
equations (once per endpoint) with **opposite signs** (reference direction
points out of one endpoint, into the other) — it cancels. What survives:
external legs attached to $A$ (i.e. $S_A$, by definition) and the single
cut edge $e$, whose *other* endpoint $v$ lies outside $A$, so only one of
its two vertex-contributions is present in the sum — it does not cancel.
Hence

$$\sum_{i\in S_A}(\text{signed})\,p_i \;+\;(\text{signed})\,q_e = 0
\quad\Longrightarrow\quad q_e = \pm\sum_{i\in S_A}p_i.$$

**Step 3 (generality).** Nothing above used the specific shape of the
tree — only that it *is* a tree (unique paths, $V{-}1$ edges) and that
conservation holds at each vertex. So this holds for every internal line
of every tree diagram contributing to the process, which is exactly the
claim used in [Topic 1](topic1-peaks-and-O2n.md)'s subset-counting argument. $\blacksquare$

### Diagram

<img src="../figures/multiperipheral.svg" alt="multiperipheral tree with q1, q2 momentum-cut demonstration" width="430">

(source: `../figures/multiperipheral.tex`, `tikz-feynman`, standard
straight-fermion/wavy-photon drawing convention.) This is the same
topology as the lecture's own ASCII diagram
(p1→A→p3, p2→B→p5, A—q1—C, B—q2—C, C→p4), redrawn with the two internal
lines colour-coded and the caption stating each cut's resulting partition
explicitly.

### Verification against our own diagram (algebraic)

Vertex conservation for the three vertices, using the sign convention
"incoming momenta positive":

- Vertex $A$ (p1 in, p3 out, q1 out): $q_1 = p_1 - p_3$
- Vertex $B$ (p2 in, p5 out, q2 out): $q_2 = p_2 - p_5$
- Vertex $C$ (q1 in, q2 in, p4 out): $q_1 + q_2 = p_4$

Checked symbolically (`sympy`) that vertex $C$'s equation is *not*
independent — it follows automatically from vertex $A$'s and vertex $B$'s
equations plus overall conservation $p_1+p_2=p_3+p_4+p_5$:

```python
import sympy as sp
p1,p2,p3,p4,p5 = sp.symbols('p1 p2 p3 p4 p5')
q1 = p1 - p3
q2 = p2 - p5
p4_expr = sp.solve(sp.Eq(p1+p2, p3+p4+p5), p4)[0]
sp.simplify((q1+q2) - p4_expr)   # -> 0
```

This is itself a direct illustration of Step 2's telescoping: summing the
conservation laws at $A$ and $B$ already forces $C$'s conservation law to
hold, exactly as the proof says redundant internal cancellation must
happen for any sub-collection of vertices.

### Verification against our own diagram (graph-theoretic)

Built the actual graph (vertices $A,B,C$; external legs as pendant nodes)
and cut each internal edge in turn with `networkx.connected_components`:

```
Cutting q1 (edge A–C): partitions external legs into {p1,p3} | {p2,p4,p5}
Cutting q2 (edge B–C): partitions external legs into {p2,p5} | {p1,p3,p4}
```

Matches both the algebra above and `pickin.c:16`'s own variable
definitions exactly: `t1=(p1-p3)^2` (the invariant built from the
$\{p_1,p_3\}$ side of the $q_1$ cut) and `t2=(p2-p5)^2` (from the
$\{p_2,p_5\}$ side of the $q_2$ cut) — i.e. `pickin.c` is choosing, for
each propagator, the *smaller/simpler* side of its cut as the invariant's
defining subset, which is always a valid and equivalent choice per the
proof ($q_e=+\sum_{S_A}p_i=-\sum_{S_B}p_i$, both sides give the same
$q_e^2$).
