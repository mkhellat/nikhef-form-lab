#!/usr/bin/env python3
"""
Ward-identity check on the SAME 16-term squared-trace machinery as
mumu_naive_trace.py / mumu_trace_lib.py.

Test (non-trivial):
  Tr[(p6-slash + m) N(e1 -> q1) (p7-slash - m) Nbar(physical e1)] == 0

Substitute e1 -> q1 only on the N side. Replacing e1 on both N and Nbar
would give 0*0=0 and prove nothing.

Diagram-specific exact replacements (from momentum routing):
  diagram 1 (N1, N2): q1 = p6 - k1
  diagram 2 (N3, N4): q1 = p7 + k2

Also substitutes D1 = k1.k1 - m^2, D2 = k2.k2 - m^2, and the on-shell
conditions p6.p6 = p7.p7 = m^2. Omitting the on-shell conditions makes
this check fail (nonzero remainder) — that was caught during development
and is documented in the report appendix.
"""
import sympy as sp

from mumu_trace_lib import (
    D1, D2, NBAR_PIECES, dot, m, on_shell_muon_subs, sandwich_trace,
)

# N with e1 -> q1 expanded into concrete chains (q1 is a sum, not a token).
N_PIECES_GAUGE = {
    'N1': [  # (p6 - k1) k1 e2 / D1
        (sp.Integer(1) / D1, ['p6', 'k1', 'e2']),
        (sp.Integer(-1) / D1, ['k1', 'k1', 'e2']),
    ],
    'N2': [  # m (p6 - k1) e2 / D1
        (m / D1, ['p6', 'e2']),
        (-m / D1, ['k1', 'e2']),
    ],
    'N3': [  # e2 k2 (p7 + k2) / D2
        (sp.Integer(1) / D2, ['e2', 'k2', 'p7']),
        (sp.Integer(1) / D2, ['e2', 'k2', 'k2']),
    ],
    'N4': [  # m e2 (p7 + k2) / D2
        (m / D2, ['e2', 'p7']),
        (m / D2, ['e2', 'k2']),
    ],
}


def main():
    S_total = sp.Integer(0)
    for i, chain_terms in N_PIECES_GAUGE.items():
        for j, (cj, Dj, mj) in NBAR_PIECES.items():
            for coeff_i, chain_i in chain_terms:
                # sandwich_trace does not attach Dj; mj is the mass factor of Nbar_j.
                # Denominator of Nbar piece is Dj; coeff_i already carries 1/Di.
                T = sandwich_trace(chain_i, sp.Integer(1), cj, mj) / Dj
                S_total += coeff_i * T

    S_total = sp.expand(S_total)
    k1k1, k2k2 = sp.symbols('k1.k1 k2.k2')
    S_total = sp.together(S_total)
    S_total = S_total.subs({D1: k1k1 - m**2, D2: k2k2 - m**2})
    S_total = S_total.subs(on_shell_muon_subs())
    S_total = sp.simplify(S_total)

    print('Tr[(p6sl+m) N(e1->q1) (p7sl-m) Nbar(physical e1)] =')
    print()
    print(f'  {S_total}')
    print()
    if S_total == 0:
        print('VANISHES -- confirms M(e1->q1)=0 propagates correctly through')
        print('the full 16-term squared-trace machinery used in Section 3.')
    else:
        print('DOES NOT VANISH -- discrepancy between the amplitude-level gauge')
        print('check and the squared-trace expression; needs investigation')
        print('before Section 3\'s result can be trusted.')
        raise SystemExit(1)


if __name__ == '__main__':
    main()
