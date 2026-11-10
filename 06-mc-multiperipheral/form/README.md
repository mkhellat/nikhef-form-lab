# FORM checks (muon-pair cancellation report)

Architecture choice for Section 5: **FORM**, not Python/C/Rust.
Dirac traces are what FORM is for; C owns numeric digit-loss
(Sections 6–7); the (3.5)/$M_{ij}$ *numeric* identities are checked in
Python (`scripts/mumu_stable_35_check.py`, `scripts/mumu_mij_check.py`)
because they need Levi-Civita bilinears against the already-generated
naive poly. FORM owns the symbolic App-C reductions and the $M_{ij}$
Levi-Civita↔poly transcription.

| File | Role | Pass criterion |
|---|---|---|
| `mumu_ward.frm` | Squared Ward identity (`e1→q1` on `N` only) | prints `Sg = 0;` |
| `mumu_naive.frm` | Naive spin-summed `\|M\|^2` via `Trace4` | runs; FORM baseline |
| `mumu_schouten_c1.frm` | Appendix-C toolkit: εε norm, biline, Schouten rewrite | `NnormCheck=Biline=SchCheck=0` |
| `mumu_schouten_c2_c14.frm` | App-C reductions (C.4)/(C.8)/(C.9)/(C.13)/(C.14) | `DiffC4…DiffC14=0` |
| `mumu_mij.frm` | (3.11) $M_{ij}$ Levi-Civita form vs closed poly | `DiffM11…DiffM22=0` |
| `Makefile` | `make check` = Ward + App-C + $M_{ij}$ | exit 0 |

```bash
make check          # Ward + schouten + appC + mij
make schouten       # Appendix-C toolkit only
make appC           # (C.4)/(C.8)/(C.9)/(C.13)/(C.14) only
make mij            # (3.11) M_ij transcription only
make naive          # optional inspection of S
```

## Appendix-C toolkit notes

- FORM's Dimension-4 `e_*e_` normalization is `+24`; the Minkowski
  `(+---)` helper in Python uses an overall minus relative to FORM for
  the two-epsilon biline (`+2*(…)` in FORM vs `-2*(…)` in Python).
- The Schouten (C1) gate installs a rewrite procedure and checks it
  applies (`SchCheck=0`).
- `mumu_schouten_c2_c14.frm` contracts the multi-epsilon chains twice
  (Dimension 4) and checks the load-bearing reductions for $M_{12}$/$M_{21}$:
  - (C.4)/(C.7), (C.8), (C.9) with paper coefficients once `Q·1·2·R` is
    written `e_(q1,e1,Q,m)*e_(q2,e2,R,m)` (the `…,m,R` order flips (C.4)
    to −½).
  - (C.13) needs the extra `¼(ε^{q1 e1 αβ}ε_{Q q2 αβ})²` term; FORM locks
    the Tr[1·1] prefactor at +½ (paper scan −½).
  - (C.14) locks with the Tr[1·1] brace sign flipped relative to the
    paper scan.

## $M_{ij}$ transcription notes

`mumu_mij.frm` writes each $M_{ij}$ as a Levi-Civita `*lc` expression
and as a closed `*poly` (matching `eval_mij.h` / the report). Locked map:

- `Tr` matches Python/C; FORM biline $= -$ Python biline; `eps3sq` $=$
  $-$ FORM three-index `e_*e_`; Python `chain3` $= -$ FORM
  `e_(p1,q1,Q,m)*e_(p2,q2,m,Q)`; four-index `e_*e_` $=$ det Gram.

## What is / is not claimed here

- **Claimed:** the squared Ward premise (`Sg = 0`) on the same
  squared-trace object Section 3 studied in Python.
- **Claimed:** Appendix-C Levi-Civita normalization + biline structure +
  Schouten rewrite procedure (`make schouten`).
- **Claimed:** FORM verification of App-C (C.4)/(C.8)/(C.9)/(C.13)/(C.14)
  (`make appC`).
- **Claimed:** FORM $M_{ij}$ Levi-Civita ↔ poly transcription
  (`make mij`).
- **Claimed elsewhere (Python/C):** naive ≡ paper (3.5); $\sum M_{ij}/B^2$ ≡
  $L\otimes L$(naive). See report appendix “Verification of (3.5) and
  $M_{ij}$”.
