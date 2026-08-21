#!/usr/bin/env python3
"""
Full symbolic trace computation of the "naive / standard fashion" squared
amplitude for gamma*gamma* -> mu- mu+ (report/mumu_cancellation.tex, Section
3), matching twophoton.pdf's warning about evaluating the matrix element
"completely in terms of 4-vector products."

Strategy: rather than trust hand contraction of an 8-gamma-matrix trace,
represent each of N's four pieces (N1..N4, see the report) as an explicit
product of gamma "slots," expand the full trace
    Tr[(p6-slash + m) N (p7-slash - m) Nbar]
using Clifford-algebra contraction (Wick's trace theorem: iteratively
contract the first gamma matrix in a chain against every other one, with a
sign for how many gammas were skipped, replacing each contracted pair by
the metric dot product), and collect the result as a sum of scalar terms
built purely from external dot products (p6.p7, e1.e2, k1.p6, etc.) --
i.e. literally the "4-vector products" form the paper describes.

We do NOT use gamma^5 / Levi-Civita anywhere: this reaction (unlike the
pi0 case) has no axial/pseudoscalar coupling, so no epsilon tensors should
appear here -- confirmed by construction (see Section 3 discussion).
"""
import itertools
import sympy as sp

# ---------------------------------------------------------------------
# Symbols: dot products of the 6 vectors that appear (e1, e2, k1, k2, p6, p7)
# plus the mass m and propagator denominators D1, D2 (kept symbolic; not
# needed for the trace algebra itself, only reinserted at the end).
# ---------------------------------------------------------------------
vectors = ['e1', 'e2', 'k1', 'k2', 'p6', 'p7']
dots = {}
for a, b in itertools.combinations_with_replacement(vectors, 2):
    dots[(a, b)] = sp.Symbol(f'{a}.{b}')
    dots[(b, a)] = dots[(a, b)]


def dot(a, b):
    return dots[(a, b)]


def trace_chain(chain):
    """
    Tr[slash(v1) slash(v2) ... slash(vn)] via the standard recursive
    contraction formula (Wick's theorem for gamma matrices):
      Tr[g1 g2 ... gn] = sum_{k=2}^n (-1)^k (g1.gk) Tr[g2...gk-1 gk+1...gn]
    with Tr[1] = 4 (we track the overall factor of 4 at the very end) and
    any trace of an odd number of gammas = 0.
    """
    n = len(chain)
    if n == 0:
        return sp.Integer(1)          # bare Tr[1] = 4 handled by caller
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
    """Tr[slash(v1)...slash(vn)] including the overall factor of 4."""
    return sp.expand(4 * trace_chain(chain))


# ---------------------------------------------------------------------
# N's four pieces, each a chain of vector labels standing for slashed
# vectors, PLUS a scalar prefactor (1/D1, m/D1, etc, tracked separately).
# We work with "reduced" chains -- the propagator's mass term is a scalar
# m multiplying a 1-gamma or bare chain, not a slashed m.
# ---------------------------------------------------------------------
# N1 = e1-slash k1-slash e2-slash / D1         (chain, coefficient m^0)
# N2 = m * e1-slash e2-slash / D1              (chain has 2 gammas, coeff m)
# N3 = e2-slash k2-slash e1-slash / D2
# N4 = m * e2-slash e1-slash / D2
N_pieces = {
    'N1': (['e1', 'k1', 'e2'], sp.Symbol('D1')**-1, 1),
    'N2': (['e1', 'e2'],       sp.Symbol('D1')**-1, sp.Symbol('m')),
    'N3': (['e2', 'k2', 'e1'], sp.Symbol('D2')**-1, 1),
    'N4': (['e2', 'e1'],       sp.Symbol('D2')**-1, sp.Symbol('m')),
}

# Nbar reverses gamma order term by term (real e1, e2 assumed):
Nbar_pieces = {
    'N1': (['e2', 'k1', 'e1'], sp.Symbol('D1')**-1, 1),
    'N2': (['e2', 'e1'],       sp.Symbol('D1')**-1, sp.Symbol('m')),
    'N3': (['e1', 'k2', 'e2'], sp.Symbol('D2')**-1, 1),
    'N4': (['e1', 'e2'],       sp.Symbol('D2')**-1, sp.Symbol('m')),
}

m = sp.Symbol('m')


def sandwich_trace(chain_i, mass_i, chain_j, mass_j):
    """
    Computes Tr[(p6-slash + m) * (mass_i * chain_i) * (p7-slash - m) *
    (mass_j * chain_j)], expanding (p6-slash + m) and (p7-slash - m) into
    their 2x2 = 4 sub-terms (p6-slash * ... * p7-slash, p6-slash * ... *
    (-m), m * ... * p7-slash, m * ... * (-m)), each a pure gamma-chain
    trace (or, for the two "one bare m" terms, a trace of chain_i+chain_j
    alone with a scalar m prefactor -- these vanish unless the combined
    chain has even length, i.e. an ODD number of gammas from chain_i+chain_j
    survives only if... handled generically by trace_chain's odd-length
    zero rule).
    """
    total = sp.Integer(0)
    # (p6-slash)(chain_i)(p7-slash)(chain_j)
    total += full_trace(['p6'] + chain_i + ['p7'] + chain_j)
    # (p6-slash)(chain_i)(-m)(chain_j)
    total += (-m) * full_trace(['p6'] + chain_i + chain_j)
    # (m)(chain_i)(p7-slash)(chain_j)
    total += m * full_trace(chain_i + ['p7'] + chain_j)
    # (m)(chain_i)(-m)(chain_j)
    total += (-m * m) * full_trace(chain_i + chain_j)
    return sp.expand(mass_i * mass_j * total)


results = {}
S_total = sp.Integer(0)
for i, (ci, Di, mi) in N_pieces.items():
    for j, (cj, Dj, mj) in Nbar_pieces.items():
        Tij = sandwich_trace(ci, mi, cj, mj)
        Tij_full = sp.expand(Tij / (sp.Symbol('D1') if 'N1' in i or 'N2' in i else sp.Symbol('D2'))
                              / (sp.Symbol('D1') if 'N1' in j or 'N2' in j else sp.Symbol('D2')))
        # On-shell muon conditions: p6, p7 are the physical outgoing mu-,
        # mu+ momenta, so p6.p6 = p7.p7 = m^2 identically -- these must
        # NOT be left as free symbols (caught via mumu_gauge_check_squared.py,
        # which found a spurious nonzero Ward-identity result before this
        # substitution was added).
        Tij_full = Tij_full.subs({dot('p6', 'p6'): m**2, dot('p7', 'p7'): m**2})
        results[(i, j)] = Tij_full
        S_total += Tij_full

S_total = sp.expand(S_total)

print("=" * 70)
print("Naive squared-amplitude trace S = sum_{i,j} T_ij  (16 terms)")
print("=" * 70)
for (i, j), val in results.items():
    nonzero = "" if val != 0 else "   [vanishes]"
    print(f"T[{i},{j}] = {val}{nonzero}")

print()
print("=" * 70)
print("Sum of all 16 terms, fully expanded:")
print("=" * 70)
print(S_total)

print()
print(f"Number of distinct monomial terms in the fully expanded sum: "
      f"{len(sp.Add.make_args(S_total))}")
