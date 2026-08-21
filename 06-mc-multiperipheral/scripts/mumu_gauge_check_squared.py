#!/usr/bin/env python3
"""
Stronger Ward-identity cross-check, at the level of the full 16-term
squared trace (scripts/mumu_naive_trace.py), not just the amplitude N
(scripts/mumu_gauge_check.py). This is the exact expression that goes
into report/mumu_cancellation.tex Section 3 and later the numerical
demonstration -- so it is the one that actually needs validating, not
just the simpler amplitude-level identity.

The physically correct test is NOT "substitute e1->q1 everywhere and see
if the 98-term sum vanishes" (that would test M(e1->q1) * Mbar(e1->q1) =
0 * 0 = 0 trivially, proving nothing new). The correct test is:
  Tr[ (p6sl+m) N(e1->q1) (p7sl-m) Nbar ] = 0
i.e. substitute e1->q1 ONLY on the N side (bra), leaving Nbar (the
physical, real-e1 conjugate amplitude) untouched -- this is the genuine
statement M(e1->q1)=0 implies M(e1->q1)*Mbar(physical) = 0 for ANY fixed
Mbar, which is a real, nontrivial constraint on the trace structure.

We reuse the same chain/trace machinery as mumu_naive_trace.py, but build
N's chains with e1 replaced by its two diagram-specific exact expressions
(q1 = p6-k1 for the N1/N2 chains that come from diagram 1's e1-slash...,
q1 = p7+k2 for the N3/N4 chains from diagram 2's ...e1-slash), while Nbar
keeps generic e1.
"""
import itertools
import sympy as sp

vectors = ['e1', 'e2', 'k1', 'k2', 'p6', 'p7']
dots = {}
for a, b in itertools.combinations_with_replacement(vectors, 2):
    dots[(a, b)] = sp.Symbol(f'{a}.{b}')
    dots[(b, a)] = dots[(a, b)]


def dot(a, b):
    return dots[(a, b)]


def trace_chain(chain):
    n = len(chain)
    if n == 0:
        return sp.Integer(1)
    if n % 2 == 1:
        return sp.Integer(0)
    if n == 2:
        return dot(chain[0], chain[1])
    total = sp.Integer(0)
    v1 = chain[0]
    for k in range(1, n):
        sign = sp.Integer(-1) ** (k - 1)
        remaining = chain[1:k] + chain[k + 1:]
        total += sign * dot(v1, chain[k]) * trace_chain(remaining)
    return sp.expand(total)


def full_trace(chain):
    return sp.expand(4 * trace_chain(chain))


m = sp.Symbol('m')
D1, D2 = sp.symbols('D1 D2')


def sandwich_trace(chain_i, chain_j):
    total = sp.Integer(0)
    total += full_trace(['p6'] + chain_i + ['p7'] + chain_j)
    total += (-m) * full_trace(['p6'] + chain_i + chain_j)
    total += m * full_trace(chain_i + ['p7'] + chain_j)
    total += (-m * m) * full_trace(chain_i + chain_j)
    return sp.expand(total)


# N's chains with e1 -> q1 substituted (diagram-specific, exact momentum
# relations, matching mumu_gauge_check.py's verified-by-hand reduction):
#   diagram 1 (N1, N2): e1 -> q1 = p6 - k1
#   diagram 2 (N3, N4): e1 -> q1 = p7 + k2
# Represented here as explicit CHAIN EXPANSIONS (each a list of
# (coefficient, chain) pairs) since q1 is now a sum of two vectors, not a
# single token -- e.g. N1 = e1sl k1sl e2sl / D1 becomes
#   (p6sl - k1sl) k1sl e2sl / D1 = p6sl k1sl e2sl / D1 - k1sl k1sl e2sl / D1
N_pieces_gauge = {
    # N1 = e1sl k1sl e2sl / D1  ->  (p6sl-k1sl) k1sl e2sl / D1
    'N1': [(sp.Integer(1) / D1, ['p6', 'k1', 'e2']),
           (sp.Integer(-1) / D1, ['k1', 'k1', 'e2'])],
    # N2 = m e1sl e2sl / D1  ->  m (p6sl-k1sl) e2sl / D1
    'N2': [(m / D1, ['p6', 'e2']),
           (-m / D1, ['k1', 'e2'])],
    # N3 = e2sl k2sl e1sl / D2  ->  e2sl k2sl (p7sl+k2sl) / D2
    'N3': [(sp.Integer(1) / D2, ['e2', 'k2', 'p7']),
           (sp.Integer(1) / D2, ['e2', 'k2', 'k2'])],
    # N4 = m e2sl e1sl / D2  ->  m e2sl (p7sl+k2sl) / D2
    'N4': [(m / D2, ['e2', 'p7']),
           (m / D2, ['e2', 'k2'])],
}

# Nbar keeps the PHYSICAL (real e1) form, unchanged from mumu_naive_trace.py:
Nbar_pieces = {
    'N1': (['e2', 'k1', 'e1'], sp.Symbol('D1')**-1, 1),
    'N2': (['e2', 'e1'],       sp.Symbol('D1')**-1, m),
    'N3': (['e1', 'k2', 'e2'], sp.Symbol('D2')**-1, 1),
    'N4': (['e1', 'e2'],       sp.Symbol('D2')**-1, m),
}

S_total = sp.Integer(0)
for i, chain_terms in N_pieces_gauge.items():
    for j, (cj, Dj, mj) in Nbar_pieces.items():
        for coeff_i, chain_i in chain_terms:
            Tij = sandwich_trace(chain_i, cj) * mj
            S_total += coeff_i * Tij

S_total = sp.expand(S_total)

# Substitute D1 = k1.k1 - m^2, D2 = k2.k2 - m^2, AND the on-shell muon
# conditions p6.p6 = p7.p7 = m^2 (these were left as free symbols above,
# which is wrong -- p6, p7 are the physical outgoing muon momenta).
k1k1, k2k2 = sp.symbols('k1.k1 k2.k2')
S_total = sp.together(S_total)
S_total = S_total.subs({sp.Symbol('D1'): k1k1 - m**2, sp.Symbol('D2'): k2k2 - m**2})
S_total = S_total.subs({dots[('p6', 'p6')]: m**2, dots[('p7', 'p7')]: m**2})
S_total = sp.simplify(S_total)

print("Tr[(p6sl+m) N(e1->q1) (p7sl-m) Nbar(physical e1)] =")
print()
print(f"  {S_total}")
print()
if S_total == 0:
    print("VANISHES -- confirms M(e1->q1)=0 propagates correctly through")
    print("the full 16-term squared-trace machinery used in Section 3.")
else:
    print("DOES NOT VANISH -- discrepancy between the amplitude-level gauge")
    print("check and the squared-trace expression; needs investigation")
    print("before Section 3's 98-term result can be trusted.")
