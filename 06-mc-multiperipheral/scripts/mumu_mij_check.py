#!/usr/bin/env python3
"""Check Vermaseren (3.11) M11 against the W1⊗W1-dressed naive poly.

M11 is the only M_ij with no Levi-Civita tensors (pure QED block). With
Q = p6−p7 and B = (q1·q2)^2−(Q·q1)^2, the identity checked here is

  M11 / B^2  ≡  S_naive contracted with L_W1(q1) ⊗ L_W1(q2)

where L_W1(q) = 2 q² (g − qq/q²) is the W1 piece of the point-particle
tensor (3.8). (The (q1² q2²)^2 in the (3.10) denominator is already
accounted for by writing the comparison as M11/B^2 ≡ S_W11.)

Pass: relative error < 1e-5 on N random multiperipheral samples
(finite-difference Hessian extraction of the W1 contraction).

M12, M21, M22 (ε structures + beam momenta) are not claimed here.

Run from scripts/:  python3 mumu_mij_check.py
"""
from __future__ import annotations

import sys

import numpy as np
import sympy as sp

from mumu_trace_lib import compute_naive_Tij


def mdot(a: np.ndarray, b: np.ndarray) -> float:
    return float(a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3])


def M11(q1, q2, Q, m) -> float:
    q1q1, q2q2, q1q2 = mdot(q1, q1), mdot(q2, q2), mdot(q1, q2)
    Qq1, QQ = mdot(Q, q1), mdot(Q, Q)
    B = q1q2**2 - Qq1**2
    brace = B * (
        Qq1**2 - q1q1 * q2q2 + q1q2**2 - QQ * (q1q1 + q2q2 + 2 * m * m)
    ) - 2 * Qq1**2 * (q1q1 + 2 * m * m) * (q2q2 + 2 * m * m)
    return 64.0 * q1q1 * q2q2 * brace


def make_eval_S():
    _, S_sym = compute_naive_Tij(apply_onshell=True)
    syms = sorted(S_sym.free_symbols, key=str)
    S_lam = sp.lambdify(syms, S_sym, "numpy")

    def eval_S(e1, e2, q1, q2, p6, p7, m):
        k1 = p6 - q1
        k2 = p6 - q2
        D1 = mdot(k1, k1) - m * m
        D2 = mdot(k2, k2) - m * m
        vecs = {"e1": e1, "e2": e2, "k1": k1, "k2": k2, "p6": p6, "p7": p7}
        vals = []
        for s in syms:
            n = str(s)
            if n == "D1":
                vals.append(D1)
            elif n == "D2":
                vals.append(D2)
            elif n == "m":
                vals.append(m)
            else:
                A, B = n.split(".")
                vals.append(mdot(vecs[A], vecs[B]))
        return float(S_lam(*vals))

    return eval_S


def contract_L_W1(q, f_of_e) -> float:
    """H_μν L^{μν} with L = 2 q² (g − qq/q²), f(e)=H_μν e^μ e^ν."""
    basis = [np.eye(4)[i] for i in range(4)]
    H = np.zeros((4, 4))
    for mu in range(4):
        for nu in range(mu, 4):
            em, en = basis[mu], basis[nu]
            if mu == nu:
                H[mu, nu] = f_of_e(em)
            else:
                H[mu, nu] = 0.5 * (f_of_e(em + en) - f_of_e(em) - f_of_e(en))
                H[nu, mu] = H[mu, nu]
    q2 = mdot(q, q)
    L_up = 2 * q2 * np.diag([1.0, -1.0, -1.0, -1.0]) - 2 * np.outer(q, q)
    return float(np.sum(H * L_up))


def sample(rng: np.random.Generator):
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
    p6 = np.array([E, p * sth * np.cos(phi), p * sth * np.sin(phi), p * cth])
    p7 = q1 + q2 - p6
    Q = p6 - p7
    return m, q1, q2, p6, p7, Q


def main() -> int:
    eval_S = make_eval_S()
    rng = np.random.default_rng(20261108)
    N = 40
    worst = 0.0
    for i in range(N):
        m, q1, q2, p6, p7, Q = sample(rng)
        B = mdot(q1, q2) ** 2 - mdot(Q, q1) ** 2
        m11 = M11(q1, q2, Q, m)
        S_W11 = contract_L_W1(
            q2,
            lambda e2: contract_L_W1(
                q1, lambda e1: eval_S(e1, e2, q1, q2, p6, p7, m)
            ),
        )
        left = m11 / B**2
        rel = abs(left - S_W11) / max(abs(S_W11), 1e-300)
        worst = max(worst, rel)
        if rel >= 1e-5:
            print(
                f"FAIL #{i}: M11/B^2={left:.12e} S_W11={S_W11:.12e} rel={rel:.3e}",
                file=sys.stderr,
            )
            return 1
    print(f"PASS: M11/B^2 ≡ W1⊗W1(naive) at {N} samples (worst rel {worst:.3e})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
