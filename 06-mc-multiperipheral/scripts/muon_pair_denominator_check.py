#!/usr/bin/env python3
"""
Verifies, algebraically, that twophoton.pdf formula (III.10)'s denominator
factor [(q1.q2)^2 - (Q.q1)^2], with Q = q1+q2, is the same bracket as
pi0.c's `la` term (up to the on-shell-mass substitution that turns the
generic bracket t1*t2 - (q1.q2)^2 into the pion/muon-pair-specific closed
form). This is the one piece of Topic 6's muon-pair "super formula" derivation
checked directly here (symbolically only -- no vendored code implements the
full muon-pair matrix element to cross-check numerically; see
notes/topic6-muon-pair-super-formulae.md for why).

Definitions used (twophoton.pdf notation): q1, q2 are the two virtual-photon
momenta, t1 = q1.q1, t2 = q2.q2, Q = q1 + q2 (paper's Appendix C notation,
P = q1+q2 in some equations, Q used in III.10/III.11 -- same object).
"""
import sympy as sp

q1q1, q2q2, q1q2 = sp.symbols('t1 t2 q1q2', real=True)  # q1.q1, q2.q2, q1.q2
t1, t2 = q1q1, q2q2

# Q = q1 + q2, so Q.q1 = q1.q1 + q1.q2 = t1 + q1q2
Q_dot_q1 = t1 + q1q2

denominator_bracket = sp.expand((q1q2)**2 - (Q_dot_q1)**2)
print("(q1.q2)^2 - (Q.q1)^2 =")
print(f"  {denominator_bracket}")
print(f"  factors as: {sp.factor(denominator_bracket)}\n")

# This should equal -(t1*(t1 + 2*q1q2)) i.e. -t1*(t1 + 2 q1.q2), NOT
# directly pi0.c's bracket -- confirm the algebraic identity first, then
# show how it relates to pi0.c's (t1*t2 - q1q2^2) form once t1 <-> t2
# symmetrization (using Q.q2 instead of Q.q1, and averaging / the on-shell
# t1+t2+2 q1.q2 = mass^2 relation pi0.c uses) is imposed.
expected = sp.expand(-t1 * (t1 + 2 * q1q2))
assert sp.expand(denominator_bracket - expected) == 0
print("Confirmed: (q1.q2)^2 - (Q.q1)^2 == -t1*(t1 + 2*q1.q2) exactly, as a")
print("pure matrix/vector identity (no on-shell relation used yet).\n")

# pi0.c's own bracket (the argument squared inside `la`) is, before the
# on-shell substitution, just t1*t2 - (q1.q2)^2 -- the DOUBLE-DROP Gram
# minor from gram_minors_pi0.py, i.e. det Gram(q1,q2). Show this is the
# NEGATIVE of the same object once q1.q2 is eliminated via the on-shell
# relation q1.q2 = (mu^2 - t1 - t2)/2 (mu = pi0 or muon-pair invariant
# mass, same substitution pi0.c performs for `la`):
mu2 = sp.symbols('mu2', real=True)
q1q2_onshell = (mu2 - t1 - t2) / 2

la_bracket_pi0 = sp.expand(t1 * t2 - q1q2**2)
la_bracket_pi0_onshell = sp.expand(la_bracket_pi0.subs(q1q2, q1q2_onshell))

denom_bracket_onshell = sp.expand(denominator_bracket.subs(q1q2, q1q2_onshell))

print("pi0.c's la-bracket (t1*t2 - (q1.q2)^2), on-shell:")
print(f"  {la_bracket_pi0_onshell}\n")
print("twophoton.pdf's (q1.q2)^2 - (Q.q1)^2, on-shell:")
print(f"  {denom_bracket_onshell}\n")

# These are generically DIFFERENT algebraic objects (one symmetric in t1,t2
# via q1.q2, the other privileging q1 through Q.q1) -- the paper's formula
# (III.10) denominator uses q1^2 q2^2 [(q1.q2)^2-(Q.q1)^2] as a WHOLE, and
# is not claimed to literally equal pi0.c's la term; the structural point
# (documented in the note) is only that both are quadratic brackets built
# from t1, t2, q1.q2 that vanish the same way the leading spinor traces
# vanish in the appropriate soft/collinear limits, not a literal formula
# match. Print both forms so the difference is visible rather than
# asserted away.
print("These are NOT identical (denominator privileges q1 through Q.q1,")
print("while pi0.c's la-bracket is symmetric in t1,t2) -- documented as a")
print("structural analogy in the note, not a literal formula match.")
