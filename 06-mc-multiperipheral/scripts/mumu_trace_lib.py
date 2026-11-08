#!/usr/bin/env python3
"""
Shared Dirac-trace machinery for the gamma*gamma* -> mu- mu+ naive
squared-amplitude calculation used by:

  - mumu_naive_trace.py          (Section 3 of report/mumu_cancellation.tex)
  - mumu_gauge_check_squared.py  (Ward identity on the same 16-term sum)

Kept in one place so the two checks cannot silently drift apart.

Conventions (must match report Sec. 2):
  N  = e1-slash (k1-slash + m)/D1 e2-slash
     + e2-slash (k2-slash + m)/D2 e1-slash
  Nbar reverses the gamma order of each piece (real e1, e2).
  Spin sum: Tr[(p6-slash + m) N (p7-slash - m) Nbar].

What this library does NOT impose:
  - momentum conservation among {e1,e2,k1,k2,p6,p7}
  - D1 = k1.k1 - m^2 or D2 = k2.k2 - m^2
  Those are left to the calling script when needed. The naive "4-vector
  products" expansion intentionally treats the external dots as
  independent symbols (apart from optional on-shell p6.p6 = p7.p7 = m^2).
"""
import itertools
import sympy as sp

VECTORS = ['e1', 'e2', 'k1', 'k2', 'p6', 'p7']

dots = {}
for a, b in itertools.combinations_with_replacement(VECTORS, 2):
    dots[(a, b)] = sp.Symbol(f'{a}.{b}')
    dots[(b, a)] = dots[(a, b)]

m = sp.Symbol('m')
D1, D2 = sp.symbols('D1 D2')


def dot(a, b):
    return dots[(a, b)]


def trace_chain(chain):
    """
    Tr[slash(v1)...slash(vn)] / 4  via recursive contraction:
      Tr[g1...gn] = sum_{k=2}^{n} (-1)^{k-1} (g1.gk) Tr[without g1,gk]
    Odd length -> 0. Empty product -> 1 (so full_trace multiplies by 4).
    """
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
    """Tr[slash(v1)...slash(vn)] including Tr[1] = 4."""
    return sp.expand(4 * trace_chain(chain))


def sandwich_trace(chain_i, mass_i, chain_j, mass_j):
    """
    Tr[(p6-slash + m) * (mass_i * chain_i) * (p7-slash - m) * (mass_j * chain_j)]
    expanded into the four (p6-slash/+m) x (p7-slash/-m) pieces. Does not
    attach propagator denominators.
    """
    total = sp.Integer(0)
    total += full_trace(['p6'] + chain_i + ['p7'] + chain_j)
    total += (-m) * full_trace(['p6'] + chain_i + chain_j)
    total += m * full_trace(chain_i + ['p7'] + chain_j)
    total += (-m * m) * full_trace(chain_i + chain_j)
    return sp.expand(mass_i * mass_j * total)


# N = N1 + N2 + N3 + N4. Each entry: (gamma-chain, denominator, mass factor).
N_PIECES = {
    'N1': (['e1', 'k1', 'e2'], D1, sp.Integer(1)),  # e1 k1 e2 / D1
    'N2': (['e1', 'e2'],       D1, m),              # m e1 e2 / D1
    'N3': (['e2', 'k2', 'e1'], D2, sp.Integer(1)),  # e2 k2 e1 / D2
    'N4': (['e2', 'e1'],       D2, m),              # m e2 e1 / D2
}

# Nbar: reverse gamma order, same denominators / mass factors.
NBAR_PIECES = {
    'N1': (['e2', 'k1', 'e1'], D1, sp.Integer(1)),
    'N2': (['e2', 'e1'],       D1, m),
    'N3': (['e1', 'k2', 'e2'], D2, sp.Integer(1)),
    'N4': (['e1', 'e2'],       D2, m),
}


def on_shell_muon_subs():
    """p6, p7 are on-shell outgoing muon momenta."""
    return {dot('p6', 'p6'): m**2, dot('p7', 'p7'): m**2}


def compute_naive_Tij(apply_onshell=True):
    """
    Return (results_dict, S_total) where results_dict[(i,j)] is T_ij
    including its 1/(Di Dj) factor.
    """
    results = {}
    S_total = sp.Integer(0)
    for i, (ci, Di, mi) in N_PIECES.items():
        for j, (cj, Dj, mj) in NBAR_PIECES.items():
            Tij = sandwich_trace(ci, mi, cj, mj) / (Di * Dj)
            Tij = sp.expand(Tij)
            if apply_onshell:
                Tij = sp.expand(Tij.subs(on_shell_muon_subs()))
            results[(i, j)] = Tij
            S_total += Tij
    return results, sp.expand(S_total)


def monomial_count(expr):
    if expr == 0:
        return 0
    return len(sp.Add.make_args(sp.expand(expr)))
