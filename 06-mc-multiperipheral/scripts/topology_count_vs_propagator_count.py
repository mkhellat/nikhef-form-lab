#!/usr/bin/env python3
"""
Verifies NOTES.md's Topic 1 correction: "number of diagrams" and "number
of potential propagator momenta" (2^(n+1)-n-3) are different
combinatorial quantities with different growth rates, not
interchangeable stand-ins for the same O(2^n) hazard count.

Two things are computed and compared:
  (1) the propagator-subset count already derived in NOTES.md,
      2^(n_ext-1) - n_ext - 1 = 2^(n+1) - n - 3
  (2) the number of distinct unrooted cubic (trivalent) tree topologies
      on n_ext labeled leaves -- the classical double-factorial count
      (2*n_ext - 5)!!, which is exactly "how many distinct tree
      Feynman diagrams exist for an n-particle final state" (ignoring
      which specific particle sits where in the internal structure,
      i.e. counting topologies, not accounting for identical particles).

(2) grows super-exponentially and overtakes (1) already by n=5-7,
confirming these should not be treated as the same order of growth.
"""
from math import factorial


def double_factorial_odd(k: int) -> int:
    result = 1
    while k > 1:
        result *= k
        k -= 2
    return result


def propagator_subset_count(n: int) -> int:
    n_ext = n + 2
    return 2 ** (n_ext - 1) - n_ext - 1


def tree_topology_count(n: int) -> int:
    n_ext = n + 2
    if n_ext < 3:
        return 1
    return double_factorial_odd(2 * n_ext - 5)


print(f"{'n':>2} {'n_ext':>6} {'propagator-subset count':>24} {'tree topology count (2n_ext-5)!!':>34} {'topology/propagator ratio':>26}")
for n in range(1, 8):
    p = propagator_subset_count(n)
    t = tree_topology_count(n)
    ratio = t / p if p else float('inf')
    print(f"{n:>2} {n+2:>6} {p:>24} {t:>34} {ratio:>26.2f}")

n7_p = propagator_subset_count(7)
n7_t = tree_topology_count(7)
assert n7_t > 500 * n7_p, "expected topology count to dwarf propagator count by n=7"
print(f"\nAt n=7: topology count is {n7_t / n7_p:.0f}x the propagator-subset count.")
print("Confirms these are genuinely different growth rates, not interchangeable.")
