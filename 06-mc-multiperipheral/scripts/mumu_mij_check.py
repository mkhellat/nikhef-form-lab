#!/usr/bin/env python3
"""Check Vermaseren (3.11) M_ij against electron-dressed naive |M|^2.

With Q = p6−p7 and B = (q1·q2)^2−(Q·q1)^2:

  M_ij / B^2  ≡  S_naive contracted with L_i(q1,p1) ⊗ L_j(q2,p2)

where the point-particle tensor (3.8) splits as
  L_W1(q) = 2 q² (g − qq/q²),
  L_W2(q,p) = 8 r_μ r_ν,  r = p − (p·q)/q² q,
and L = L_W1 + L_W2 (W1=W2=1).

Gates:
  M11  — W1⊗W1 (no ε; already stable)
  M12  — W1⊗W2
  M21  — W2⊗W1
  M22  — W2⊗W2  (parallel to (3.5) with (qi,ei)→(pi,qi))
  full — Σ M_ij vs L⊗L

Run from scripts/:  python3 mumu_mij_check.py
"""
from __future__ import annotations

import sys

import numpy as np
import sympy as sp
from numpy.linalg import det

from mumu_trace_lib import compute_naive_Tij
from mumu_stable_35_check import Tr, biline, eps_full, mdot


def eps3sq(a, b, c) -> float:
    """Gram stand-in for ε^{a b c μ} ε_{a b c μ} in (+−−−): −det Gram.

    Relative to NPB (3.11) typesetting of the three-index contractions in
    M12/M21, this helper has the opposite overall sign (verified numerically
    against W1⊗W2 / W2⊗W1). The M12/M21 formulas below absorb that once.
    """
    G = np.array(
        [
            [mdot(a, a), mdot(a, b), mdot(a, c)],
            [mdot(b, a), mdot(b, b), mdot(b, c)],
            [mdot(c, a), mdot(c, b), mdot(c, c)],
        ]
    )
    return -float(det(G))


def chain3(a, b, c, d, e, f) -> float:
    """ε^{a b c μ} ε_{d e μ f} via 3×3 Gram (sign as in (3.5) Q·1·2·Q)."""
    M = np.array(
        [
            [mdot(a, d), mdot(a, e), mdot(a, f)],
            [mdot(b, d), mdot(b, e), mdot(b, f)],
            [mdot(c, d), mdot(c, e), mdot(c, f)],
        ]
    )
    return float(det(M))


def M11(q1, q2, Q, m) -> float:
    q1q1, q2q2, q1q2 = mdot(q1, q1), mdot(q2, q2), mdot(q1, q2)
    Qq1, QQ = mdot(Q, q1), mdot(Q, Q)
    B = q1q2**2 - Qq1**2
    brace = B * (
        Qq1**2 - q1q1 * q2q2 + q1q2**2 - QQ * (q1q1 + q2q2 + 2 * m * m)
    ) - 2 * Qq1**2 * (q1q1 + 2 * m * m) * (q2q2 + 2 * m * m)
    return 64.0 * q1q1 * q2q2 * brace


def M12(q1, q2, Q, p2, m) -> float:
    """NPB (3.11) M12; three-index εε via +eps3sq (see eps3sq docstring)."""
    q1q1, q1q2 = mdot(q1, q1), mdot(q1, q2)
    Qq1 = mdot(Q, q1)
    B = q1q2**2 - Qq1**2
    # NPB writes B*(-εε(p2,q2,q1)-εε(p2,q2,Q)); with our eps3sq that is +eps3sq.
    t1 = eps3sq(p2, q2, q1) + eps3sq(p2, q2, Q)
    e_p2q2 = biline(p2, q2, p2, q2)
    e_mix = biline(p2, q2, q1, Q)
    t2 = Qq1**2 * e_p2q2 + 0.5 * e_mix**2
    return 128.0 * q1q1 * (B * t1 - (q1q1 + 2 * m * m) * t2)


def M21(q1, q2, Q, p1, m) -> float:
    """NPB (3.11) M21; same eps3sq convention as M12."""
    q2q2, q1q2 = mdot(q2, q2), mdot(q1, q2)
    Qq1 = mdot(Q, q1)
    B = q1q2**2 - Qq1**2
    t1 = eps3sq(p1, q1, q2) + eps3sq(p1, q1, Q)
    e_p1q1 = biline(p1, q1, p1, q1)
    e_mix = biline(p1, q1, q2, Q)
    t2 = Qq1**2 * e_p1q1 + 0.5 * e_mix**2
    return 128.0 * q2q2 * (B * t1 - (q2q2 + 2 * m * m) * t2)


def M22(q1, q2, Q, p1, p2) -> float:
    """NPB (3.11) M22 = 512 × brace; bilinears are ε(p2,q2;q1,Q) etc."""
    q1q2 = mdot(q1, q2)
    Qq1 = mdot(Q, q1)
    B = q1q2**2 - Qq1**2
    Tr11 = Tr(p1, q1, p1, q1)
    Tr22 = Tr(p2, q2, p2, q2)
    Tr12 = Tr(p1, q1, p2, q2)
    EpsFull = eps_full(p1, q1, p2, q2)
    # ε^{p1 q1 Q μ} ε_{p2 q2 Q μ}
    mid_eps = chain3(p1, q1, Q, p2, q2, Q)
    mid = mid_eps - 0.5 * Tr12 * (mdot(Q, Q) + q1q2)
    # paper: (ε^{p2 q2} ε_{q1 Q})^2 and (ε^{p1 q1} ε_{q2 Q})^2 — not vs p_i
    bil_2q1Q = biline(p2, q2, q1, Q)
    bil_1q2Q = biline(p1, q1, q2, Q)
    brace = (
        B * (0.25 * Tr12**2 - EpsFull)
        - mid**2
        - 0.125 * Tr11 * bil_2q1Q**2
        - 0.125 * Tr22 * bil_1q2Q**2
        - 0.25 * Qq1**2 * Tr11 * Tr22
    )
    return 512.0 * brace


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


def L_up_matrix(q, p, which: str) -> np.ndarray:
    """Contravariant L^{μν}: 'W1', 'W2', or 'full'."""
    q2 = mdot(q, q)
    eta = np.diag([1.0, -1.0, -1.0, -1.0])
    Lw1 = 2 * q2 * eta - 2 * np.outer(q, q)
    r = p - q * (mdot(p, q) / q2)
    Lw2 = 8.0 * np.outer(r, r)
    if which == "W1":
        return Lw1
    if which == "W2":
        return Lw2
    return Lw1 + Lw2


def hessian(f_of_e) -> np.ndarray:
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
    return H


def contract(q, p, which: str, f_of_e) -> float:
    H = hessian(f_of_e)
    return float(np.sum(H * L_up_matrix(q, p, which)))


def make_massless_p(q, pT: float, rng: np.random.Generator) -> np.ndarray:
    """Massless p with 2 p·q = q² (outgoing electron after photon emission)."""
    q0, qz = q[0], q[3]
    t = mdot(q, q)
    # p = (E, pT*cosφ, pT*sinφ, pz), E² - pT² - pz² = 0, E q0 - pz qz = t/2
    phi = rng.uniform(0, 2 * np.pi)
    a = (qz / q0) ** 2 - 1.0
    b = (t * qz) / (q0**2)
    c = (t / (2 * q0)) ** 2 - pT**2
    disc = b * b - 4 * a * c
    tries = 0
    while disc < 0 and tries < 8:
        pT *= 0.5
        c = (t / (2 * q0)) ** 2 - pT**2
        disc = b * b - 4 * a * c
        tries += 1
    if disc < 0:
        raise RuntimeError("no massless p for this q,pT")
    pz = (-b + np.sqrt(disc)) / (2 * a)
    E = (t / 2 + pz * qz) / q0
    return np.array([E, pT * np.cos(phi), pT * np.sin(phi), pz])


def sample(rng: np.random.Generator):
    m = 0.105658
    E = 10.0
    epsv = 10 ** rng.uniform(-2.5, -1.0)
    qz = np.sqrt(E * E + epsv)
    pmu = np.sqrt(E * E - m * m)
    q1 = np.array([E, 0.0, 0.0, qz])
    q2 = np.array([E, 0.0, 0.0, -qz])
    cth = rng.uniform(-0.9, 0.9)
    sth = np.sqrt(max(0.0, 1.0 - cth * cth))
    phi = rng.uniform(0, 2 * np.pi)
    p6 = np.array(
        [E, pmu * sth * np.cos(phi), pmu * sth * np.sin(phi), pmu * cth]
    )
    p7 = q1 + q2 - p6
    Q = p6 - p7
    p1 = make_massless_p(q1, 0.4 + 0.4 * rng.random(), rng)
    p2 = make_massless_p(q2, 0.4 + 0.4 * rng.random(), rng)
    return m, q1, q2, p6, p7, Q, p1, p2


def dressed(eval_S, q1, q2, p6, p7, p1, p2, m, w1: str, w2: str) -> float:
    return contract(
        q2,
        p2,
        w2,
        lambda e2: contract(
            q1, p1, w1, lambda e1: eval_S(e1, e2, q1, q2, p6, p7, m)
        ),
    )


def relerr(a, b) -> float:
    return abs(a - b) / max(abs(b), 1e-300)


def main() -> int:
    eval_S = make_eval_S()
    rng = np.random.default_rng(20261108)
    N = 30
    worst = {k: 0.0 for k in ("M11", "M12", "M21", "M22", "full")}

    for i in range(N):
        try:
            m, q1, q2, p6, p7, Q, p1, p2 = sample(rng)
        except RuntimeError:
            continue
        # sanity: massless + 2p·q=q²
        if abs(mdot(p1, p1)) > 1e-8 or abs(mdot(p2, p2)) > 1e-8:
            continue
        if abs(2 * mdot(p1, q1) - mdot(q1, q1)) > 1e-6:
            continue
        if abs(2 * mdot(p2, q2) - mdot(q2, q2)) > 1e-6:
            continue

        B = mdot(q1, q2) ** 2 - mdot(Q, q1) ** 2
        checks = {
            "M11": (
                M11(q1, q2, Q, m) / B**2,
                dressed(eval_S, q1, q2, p6, p7, p1, p2, m, "W1", "W1"),
            ),
            "M12": (
                M12(q1, q2, Q, p2, m) / B**2,
                dressed(eval_S, q1, q2, p6, p7, p1, p2, m, "W1", "W2"),
            ),
            "M21": (
                M21(q1, q2, Q, p1, m) / B**2,
                dressed(eval_S, q1, q2, p6, p7, p1, p2, m, "W2", "W1"),
            ),
            "M22": (
                M22(q1, q2, Q, p1, p2) / B**2,
                dressed(eval_S, q1, q2, p6, p7, p1, p2, m, "W2", "W2"),
            ),
        }
        mij_sum = (
            M11(q1, q2, Q, m)
            + M12(q1, q2, Q, p2, m)
            + M21(q1, q2, Q, p1, m)
            + M22(q1, q2, Q, p1, p2)
        )
        checks["full"] = (
            mij_sum / B**2,
            dressed(eval_S, q1, q2, p6, p7, p1, p2, m, "full", "full"),
        )

        for name, (left, right) in checks.items():
            r = relerr(left, right)
            worst[name] = max(worst[name], r)
            tol = 1e-5 if name == "M11" else 1e-4
            if r >= tol:
                print(
                    f"FAIL {name}#{i}: left={left:.12e} right={right:.12e} "
                    f"rel={r:.3e}",
                    file=sys.stderr,
                )
                return 1

    parts = ", ".join(f"{k} {worst[k]:.3e}" for k in worst)
    print(f"PASS: M_ij/B^2 ≡ W_i⊗W_j(naive) at {N} samples (worst: {parts})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
