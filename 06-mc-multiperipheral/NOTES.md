# Discussion notes — Part 5 (multiperipheral kinematics)

Running log of questions/derivations as we work through uam19 `part5.pdf`
together, plus the underlying primary sources it cites (J.A.M. Vermaseren,
*Nucl. Phys.* B229 (1983) 347–371; E. Byckling & K. Kajantie, *Particle
Kinematics* (1973); N. Byers & C.N. Yang, *Rev. Mod. Phys.* 36 (1964) 595).

This file is an index. Each topic below has its own full discussion file
under `notes/` — the derivation, worked examples, code citations, and (where
relevant) the corrections made along the way, kept visible rather than
silently overwritten. `README.md` summarizes the settled conclusions once
enough ground is covered.

## Topics

- [Topic 0: what "tree diagram" and "loop" actually mean](notes/topic0-tree-vs-loop.md)
  — constraint-counting definition ($L=E-(V-1)$), three worked examples,
  two corrections made along the way.
- [Topic 1: "peaks", "types of peaks", and why propagator count grows like $O(2^n)$](notes/topic1-peaks-and-O2n.md)
  — the $2^{n+1}-n-3$ tree-restricted count, the loop/Landau-singularity
  generalization, why splitting phase space addresses the mismatch
  (grounded in Vermaseren's own two-method taxonomy), and the corrected
  Byers–Yang citation.
- [Topic 2: "the number of non-trivial integration variables is at most $3n-4$"](notes/topic2-3n-minus-4.md)
  — the $3n-4$/$3n-5$ accounting, checked directly against `pickin.c`/
  `orient.c`/`eepi.c` and against Byckling–Kajantie's Table III.1.
- [Topic 3: internal-line momentum is a signed subset-sum of external momenta](notes/topic3-subset-sum-proof.md)
  — the rigorous tree-graph proof underlying Topic 1's counting argument.
- [Topic 4: why the pseudoscalar ($\pi^0$) coupling is the cleanest worked example](notes/topic4-pseudoscalar-coupling.md)
  — why Vermaseren introduces the numerical-stability method via
  $e^-e^+\to e^-e^+\pi^0$ rather than the scalar or QED muon-pair case.
- [Topic 5: what "would grow with powers of $s$" actually refers to in formula (3.2)](notes/topic5-s-scaling-pi0.md)
  — why $s$-growth at fixed $t_1,t_2$ would violate unitarity, the
  analytic $s^2$ derivation of the dominant Gram-determinant term, and
  its numerical confirmation directly against `pickin`/`orient`/`pi0`
  (no cancellation between `part1..part4`, correcting an earlier
  assumption in Topic 4).

## Open questions not yet filed into a topic

(New questions get discussed here first, then moved into a topic file — new
or existing — once the answer settles.)
