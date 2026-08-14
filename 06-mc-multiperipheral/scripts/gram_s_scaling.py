#!/usr/bin/env python3
"""
Derives the s-scaling of the Gram determinant eps^{p1 q1 p2 q2}
eps_{p1 q1 p2 q2} (pi0.c's `levi.gram`, formula (3.2)'s first/dominant
term) at FIXED small photon virtualities t1=q1^2, t2=q2^2, as the beam
energy sqrt(s) -> infinity. This is the numerator that the uam19 part5.pdf
lecture warns about: "the first term combines with 1/(t1 t2)^2 which by
itself would cause problems with physics because it would give a
crosssection that grows with powers of s."

Method: build explicit 4-vectors with the physically correct leading-order
scaling in the near-real-photon (equivalent-photon) kinematic regime --
beam momenta p1,p2 grow like sqrt(s), while the photon momenta q1,q2 stay
O(1) (their transverse component is fixed by t1,t2, and their component
along the beam direction is suppressed by an extra power of 1/sqrt(s), so
at LEADING order in 1/sqrt(s) they can be treated as fixed, purely
transverse vectors). Then compute the Gram (dot-product) determinant of
(p1,q1,p2,q2) symbolically in sqrt(s) and read off the power.

This is a *leading-order-in-1/sqrt(s)* scaling check, not a full on-shell
kinematics generator (that is done separately, exactly, using the real
pickin/orient code -- see scripts/scan_pi0_scaling.c, whose numerical
output this script's symbolic prediction should match).
"""
import sympy as sp

sqrts, qT1, qT2 = sp.symbols('sqrts qT1 qT2', positive=True)

E = sqrts / 2
# Beam momenta: exact leading order in sqrt(s) (mass terms are O(1/sqrt(s))
# corrections, dropped here since we only want the LEADING power of sqrts).
p1 = sp.Matrix([E, 0, 0, E])
p2 = sp.Matrix([E, 0, 0, -E])

# Photon momenta at fixed t1 = q1^2 = -qT1^2, t2 = q2^2 = -qT2^2: to leading
# order in 1/sqrt(s) they are purely transverse and independent of sqrt(s)
# (their energy/longitudinal pieces are down by a relative O(1/sqrt(s)),
# which is exactly the statement "t1, t2 fixed and small while s -> infinity").
q1 = sp.Matrix([0, qT1, 0, 0])
q2 = sp.Matrix([0, 0, qT2, 0])


def dot(a, b):
    """Mostly-plus-minus (+,-,-,-) Minkowski dot product."""
    return a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]


vecs = [p1, q1, p2, q2]
gram_matrix = sp.Matrix(4, 4, lambda i, j: dot(vecs[i], vecs[j]))
gram_det = sp.simplify(gram_matrix.det())

print("Gram matrix of (p1, q1, p2, q2), leading order in 1/sqrt(s):")
sp.pprint(gram_matrix)
print()
print("Gram determinant (= eps^{p1 q1 p2 q2} eps_{p1 q1 p2 q2}, up to sign/normalization):")
print(f"  {gram_det}")

degree = sp.Poly(gram_det, sqrts).degree()
print(f"\nLeading power of sqrt(s): {degree}  =>  power of s: {degree // 2}")

t1_val, t2_val = -qT1**2, -qT2**2
gram_in_t = gram_det.subs({qT1: sp.sqrt(-t1_val), qT2: sp.sqrt(-t2_val)})
print(f"\nIn terms of t1 = -qT1^2, t2 = -qT2^2: Gram ~ (t1 * t2) * s^2")
print("  (matches sqrts^4 = s^2, with qT1^2 qT2^2 = t1 t2 up to sign)")

expected_degree = 4  # s^2 in sqrt(s)^4
assert degree == expected_degree, (
    f"expected sqrt(s)^{expected_degree} (i.e. s^2) leading behavior, got sqrt(s)^{degree}"
)
print("\nConfirmed: Gram determinant ~ s^2 * t1 * t2 at fixed small t1, t2.")
print("Combined with part1 = -64*Gram/(t1 t2)^2 from pi0.c, this gives")
print("part1 ~ s^2 / (t1 t2)  -- i.e. GROWS like s^2 at fixed t1, t2, exactly")
print("the growth the lecture note warns about for the naively-expected-small")
print("numerator. See scripts/scan_pi0_scaling.c for the numerical confirmation")
print("that this growth survives, unchanged in order, in the FULL pi0() sum")
print("part1+part2+part3+part4 -- there is no cancellation between the four")
print("Levi-Civita terms of formula (3.2) (the paper states this explicitly,")
print("p.353: 'all 4 terms are positive so no cancellations occur').")
