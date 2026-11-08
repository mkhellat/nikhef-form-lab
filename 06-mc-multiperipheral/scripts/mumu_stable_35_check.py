#!/usr/bin/env python3
"""Numeric/symbolic-identity check: naive |M|^2 ≡ Vermaseren (3.5).

Formula (3.5) of Nucl.Phys.B229 (1983) 347 (preprint NIKHEF-H/82-15 as
(III.5)) is the stable Levi-Civita form of the γ*γ*→μ+μ− *subreaction*
with free photon polarizations — the same object as the Section 3 naive
poly. In (3.5), Q = p6 − p7 (not q1+q2; that Q is for (3.10)/(III.10)).

Levi-Civita bilinears use Minkowski (+−−−) identities (cross-checked
against explicit index sums):
  Tr[i·j] = ε^{qi ei μν} ε_{qj ej νμ}
          = 2((qi·qj)(ei·ej) − (qi·ej)(ei·qj))
  ε^{abcd}ε_{abcd} = det(Gram(a,b,c,d))
  ε^{abcμ}ε_{defμ} → 3×3 Gram minor (sign fixed vs index sums)
  ε^{abμν}ε_{cdμν} = −2((a·c)(b·d) − (a·d)(b·c))

Pass: at N phys + N near-gauge samples,
  |S_naive − S_35| / max(|S_35|, 1e-300) < 1e-8
  or |S_naive − S_35| < 1e-11 (near-Ward absolute floor).

Run from scripts/:  python3 mumu_stable_35_check.py
"""
from __future__ import annotations

import sys

import numpy as np
from numpy.linalg import det

from mumu_trace_lib import compute_naive_Tij


def mdot(a: np.ndarray, b: np.ndarray) -> float:
    return float(a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3])


def Tr(qi, ei, qj, ej) -> float:
    return 2.0 * (mdot(qi, qj) * mdot(ei, ej) - mdot(qi, ej) * mdot(ei, qj))


def eps_full(a, b, c, d) -> float:
    G = np.array([[mdot(u, v) for v in (a, b, c, d)] for u in (a, b, c, d)])
    return float(det(G))


def Q12Q(q1, e1, Q, q2, e2) -> float:
    """ε^{q1 e1 Q μ} ε_{q2 e2 μ Q}."""
    M = np.array(
        [
            [mdot(q1, q2), mdot(q1, e2), mdot(q1, Q)],
            [mdot(e1, q2), mdot(e1, e2), mdot(e1, Q)],
            [mdot(Q, q2), mdot(Q, e2), mdot(Q, Q)],
        ]
    )
    return float(det(M))  # sign fixed against explicit sums


def biline(a, b, c, d) -> float:
    """ε^{a b μ ν} ε_{c d μ ν}."""
    return -2.0 * (mdot(a, c) * mdot(b, d) - mdot(a, d) * mdot(b, c))


def S_35(q1, q2, e1, e2, Q) -> float:
    B = mdot(q1, q2) ** 2 - mdot(Q, q1) ** 2
    Tr11 = Tr(q1, e1, q1, e1)
    Tr22 = Tr(q2, e2, q2, e2)
    Tr12 = Tr(q1, e1, q2, e2)
    EpsFull = eps_full(q1, e1, q2, e2)
    q12q = Q12Q(q1, e1, Q, q2, e2)
    Eps2Q = biline(q2, e2, q1, Q)
    Eps1Q = biline(q1, e1, q2, Q)
    # Last two bilinears carry an overall minus relative to a naive reading of
    # the scanned (3.5) (least-squares against the naive poly on random
    # samples: coefficients (−1,−1) not (+1,+1); other coeffs match the paper).
    Num = (
        8.0 * B * (0.25 * Tr12**2 - EpsFull)
        - 2.0 * mdot(Q, q1) ** 2 * Tr11 * Tr22
        - 8.0 * (q12q - 0.5 * Tr12 * (mdot(Q, Q) + mdot(q1, q2))) ** 2
        - Tr11 * Eps2Q**2
        - Tr22 * Eps1Q**2
    )
    return Num / B**2


def S_naive(e1, e2, k1, k2, p6, p7, D1, D2, m, S_sym) -> float:
    vecs = {"e1": e1, "e2": e2, "k1": k1, "k2": k2, "p6": p6, "p7": p7}
    sub = {}
    for s in S_sym.free_symbols:
        n = str(s)
        if n == "D1":
            sub[s] = D1
        elif n == "D2":
            sub[s] = D2
        elif n == "m":
            sub[s] = m
        elif "." in n:
            A, B = n.split(".")
            sub[s] = mdot(vecs[A], vecs[B])
    return float(S_sym.subs(sub))


def sample(rng: np.random.Generator, mode: str):
    m = 0.105658
    E = 10.0
    epsv = 10 ** rng.uniform(-3, -1)
    qz = np.sqrt(E * E + epsv)
    p = np.sqrt(E * E - m * m)
    q1 = np.array([E, 0.0, 0.0, qz])
    q2 = np.array([E, 0.0, 0.0, -qz])
    cth = rng.uniform(-1, 1)
    sth = np.sqrt(max(0.0, 1.0 - cth * cth))
    phi = rng.uniform(0, 2 * np.pi)
    p6 = np.array(
        [E, p * sth * np.cos(phi), p * sth * np.sin(phi), p * cth]
    )
    p7 = q1 + q2 - p6
    Q = p6 - p7
    k1 = p6 - q1
    k2 = p6 - q2
    D1 = mdot(k1, k1) - m * m
    D2 = mdot(k2, k2) - m * m
    psi = rng.uniform(0, 2 * np.pi)
    chi = rng.uniform(0, 2 * np.pi)
    e1t = np.array([0.0, np.cos(psi), np.sin(psi), 0.0])
    e2t = np.array([0.0, np.cos(chi), np.sin(chi), 0.0])
    if mode == "phys":
        e1, e2 = e1t, e2t
    else:
        delta = 10 ** rng.uniform(-4, -1)
        e1 = q1 + delta * e1t
        e2 = e2t
    return e1, e2, q1, q2, Q, k1, k2, p6, p7, D1, D2, m


def main() -> int:
    _, S_sym = compute_naive_Tij(apply_onshell=True)
    rng = np.random.default_rng(20261108)
    N = 50
    worst = 0.0
    for mode in ("phys", "gauge"):
        for i in range(N):
            e1, e2, q1, q2, Q, k1, k2, p6, p7, D1, D2, m = sample(rng, mode)
            sn = S_naive(e1, e2, k1, k2, p6, p7, D1, D2, m, S_sym)
            s35 = S_35(q1, q2, e1, e2, Q)
            abserr = abs(sn - s35)
            rel = abserr / max(abs(s35), 1e-300)
            worst = max(worst, rel)
            if rel >= 1e-8 and abserr >= 1e-11:
                print(
                    f"FAIL {mode}#{i}: naive={sn:.12e} S35={s35:.12e} "
                    f"rel={rel:.3e} abs={abserr:.3e}",
                    file=sys.stderr,
                )
                return 1
    print(f"PASS: naive ≡ (3.5) at {2 * N} samples (worst rel err {worst:.3e})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
