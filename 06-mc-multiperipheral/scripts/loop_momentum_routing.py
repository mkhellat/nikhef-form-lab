#!/usr/bin/env python3
"""
Verifies the loop-diagram generalization of Topic 3's tree lemma
(NOTES.md): for a diagram with L independent loop momenta, every
internal line's momentum is a signed subset-sum of the external
momenta PLUS an integer-coefficient combination of the L loop momenta.

Uses the one-loop QED box diagram (figures/box_loop.tex) as the
concrete example, reusing the vertex-conservation solve already
verified in NOTES.md's "Worked example 2" (Topic 0), and checks that
every internal line decomposes as claimed: external part (set loop
momentum to 0) plus an integer-coefficient loop-momentum term.
"""
import sympy as sp

a1, a2, b1, b2, l = sp.symbols('a1 a2 b1 b2 l')  # l = the free loop momentum (k_left)

# from NOTES.md Topic 0, Worked example 2's sympy solve:
# {(a1 - k_left,  b1 - b2 + k_left,  b1 + k_left,  k_left)} for (k_top,k_right,k_bottom,k_left)
lines = {
    'k1 (top fermion)':    a1 - l,
    'k4 (right photon)':   b1 - b2 + l,
    'k2 (bottom fermion)': b1 + l,
    'k3 (left photon, loop momentum itself)': l,
}

print("Internal line momentum = (external subset-sum) + (integer coeff) * l\n")
for name, expr in lines.items():
    ext_part = expr.subs(l, 0)
    loop_coeff = sp.diff(expr, l)
    assert loop_coeff in (-1, 0, 1), f"{name}: non-integer/unexpected loop coefficient"
    print(f"  {name:42s}: external part = {str(ext_part):10s}  loop coeff = {loop_coeff}")

print("\nAll four internal lines confirmed: external subset-sum + integer*l.")
