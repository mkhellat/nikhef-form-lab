#!/usr/bin/env python3
"""
Ward-identity / gauge-invariance check on the gamma*gamma* -> mu- mu+
amplitude N (report/mumu_cancellation.tex, Sections 2-3), done BEFORE
squaring -- the physically correct place to test gauge invariance, since
it is a statement about M itself (M(e1->q1) = 0), not about |M|^2.

N = e1-slash (k1-slash+m)/D1 e2-slash + e2-slash (k2-slash+m)/D2 e1-slash

with k1 = p6 - q1, k2 = p6 - q2, D1 = k1.k1 - m^2, D2 = k2.k2 - m^2, and
momentum conservation q1 + q2 = p6 + p7 (so q1 = p6 - k1 = p7 + k2, two
exact, independent ways to eliminate q1 -- one per diagram).

Reduction is done symbolically as noncommuting "slashed vector" tokens,
with two on-shell rewrite rules applied whenever p6 sits at the LEFT end
of a word (ubar(p6) p6-slash ... = m ubar(p6) ...) or p7 sits at the
RIGHT end (... p7-slash v(p7) = -m ... v(p7)), plus the Clifford identity
k-slash k-slash = k.k (a scalar) for a repeated adjacent token.

Diagram 1 (e1 -> q1 = p6 - k1):
  (p6sl - k1sl)(k1sl+m) e2sl / D1
  = [ p6sl k1sl e2sl + m p6sl e2sl - k1sl k1sl e2sl - m k1sl e2sl ] / D1
  = [ m k1sl e2sl + m^2 e2sl - k1.k1 e2sl - m k1sl e2sl ] / D1   (p6sl->m)
  = (m^2 - k1.k1) e2sl / D1  =  -D1 e2sl / D1  =  -e2sl

Diagram 2 (e1 -> q1 = p7 + k2):
  e2sl (k2sl+m) (p7sl+k2sl) / D2
  = [ e2sl k2sl p7sl + e2sl k2sl k2sl + m e2sl p7sl + m e2sl k2sl ] / D2
  = [ -m e2sl k2sl + k2.k2 e2sl - m^2 e2sl + m e2sl k2sl ] / D2  (p7sl->-m)
  = (k2.k2 - m^2) e2sl / D2  =  D2 e2sl / D2  =  +e2sl

Sum = -e2sl + e2sl = 0: gauge invariance holds only for the SUM of the two
diagrams, confirming the paper's statement that neither diagram is
separately gauge invariant.

This script reproduces that reduction generically (not hard-coding the
m*k1sl*e2sl cancellation, but deriving it from the same token-rewrite
rules used elsewhere in this directory's checks), as a machine-checked
version of the by-hand derivation above.
"""
import sympy as sp

dots = {}


def dot(a, b):
    key = tuple(sorted([a, b]))
    if key not in dots:
        dots[key] = sp.Symbol(f'{key[0]}.{key[1]}')
    return dots[key]


def move_p6_to_front(word):
    """One anticommutation step moving a 'p6' token toward index 0."""
    i = word.index('p6')
    if i == 0:
        return [(sp.Symbol('m'), word[1:])]
    x = word[i - 1]
    return [
        (2 * dot(x, 'p6'), word[:i - 1] + word[i + 1:]),
        (sp.Integer(-1), word[:i - 1] + ['p6', x] + word[i + 1:]),
    ]


def move_p7_to_back(word):
    """One anticommutation step moving a 'p7' token toward the last index."""
    i = word.index('p7')
    if i == len(word) - 1:
        return [(-sp.Symbol('m'), word[:-1])]
    y = word[i + 1]
    return [
        (2 * dot('p7', y), word[:i] + word[i + 2:]),
        (sp.Integer(-1), word[:i] + [y, 'p7'] + word[i + 1:]),
    ]


def collapse_repeats(word, coeff):
    """k-slash k-slash (adjacent identical tokens) -> k.k (a scalar)."""
    w = list(word)
    changed = True
    while changed:
        changed = False
        for i in range(len(w) - 1):
            if w[i] == w[i + 1]:
                coeff = coeff * dot(w[i], w[i])
                w = w[:i] + w[i + 2:]
                changed = True
                break
    return tuple(w), coeff


def reduce_word(coeff, word, max_iter=50):
    terms = [(coeff, list(word))]
    for _ in range(max_iter):
        new_terms, changed = [], False
        for c, w in terms:
            if 'p6' in w:
                changed = True
                new_terms += [(sp.expand(c * c2), w2) for c2, w2 in move_p6_to_front(w)]
            elif 'p7' in w:
                changed = True
                new_terms += [(sp.expand(c * c2), w2) for c2, w2 in move_p7_to_back(w)]
            else:
                new_terms.append((c, w))
        terms = new_terms
        if not changed:
            break
    collected = {}
    for c, w in terms:
        key, c = collapse_repeats(w, c)
        collected[key] = collected.get(key, sp.Integer(0)) + c
    return {k: sp.expand(v) for k, v in collected.items() if sp.expand(v) != 0}


D1, D2, m = sp.symbols('D1 D2 m')
total = {}


def add(coeff, word):
    for w, c in reduce_word(coeff, word).items():
        total[w] = total.get(w, sp.Integer(0)) + c


# Diagram 1: e1sl -> q1sl = p6sl - k1sl, giving (p6sl-k1sl)(k1sl+m)e2sl/D1
#   = [p6sl k1sl e2sl + m p6sl e2sl - k1sl k1sl e2sl - m k1sl e2sl] / D1
add(sp.Integer(1) / D1, ['p6', 'k1', 'e2'])
add(m / D1, ['p6', 'e2'])
add(sp.Integer(-1) / D1, ['k1', 'k1', 'e2'])
add(-m / D1, ['k1', 'e2'])

# Diagram 2: e1sl -> q1sl = p7sl + k2sl, giving e2sl(k2sl+m)(p7sl+k2sl)/D2
#   = [e2sl k2sl p7sl + e2sl k2sl k2sl + m e2sl p7sl + m e2sl k2sl] / D2
add(sp.Integer(1) / D2, ['e2', 'k2', 'p7'])
add(sp.Integer(1) / D2, ['e2', 'k2', 'k2'])
add(m / D2, ['e2', 'p7'])
add(m / D2, ['e2', 'k2'])

# Substitute the propagator denominators' actual definitions.
k1k1, k2k2 = sp.symbols('k1.k1 k2.k2')
total = {w: sp.simplify(c.subs({D1: k1k1 - m**2, D2: k2k2 - m**2}))
         for w, c in total.items()}
total = {w: c for w, c in total.items() if c != 0}

print("Terms remaining after e1 -> q1, on-shell Dirac equations, and")
print("D1 = k1.k1 - m^2, D2 = k2.k2 - m^2 substituted:")
print()
if not total:
    print("  ALL TERMS CANCEL -- gauge invariance (Ward identity) confirmed.")
    print("  Matches the by-hand derivation in this script's docstring:")
    print("  diagram 1 -> -e2-slash, diagram 2 -> +e2-slash, sum = 0.")
else:
    for w, c in total.items():
        print(f"  word={w}  coeff={c}")
    print()
    print("NOT ALL TERMS CANCELED -- discrepancy with the hand derivation,")
    print("needs further debugging before Section 3 can be trusted.")
