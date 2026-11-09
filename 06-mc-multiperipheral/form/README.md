# FORM checks (muon-pair cancellation report)

Architecture choice for Section 5: **FORM**, not Python/C/Rust.
Dirac traces are what FORM is for; C owns numeric digit-loss
(Sections 6–7); the (3.5)/$M_{ij}$ identities are checked in Python
(`scripts/mumu_stable_35_check.py`, `scripts/mumu_mij_check.py`) because
they need Levi-Civita bilinears against the already-generated naive poly.
Rust would not buy domain leverage here.

| File | Role | Pass criterion |
|---|---|---|
| `mumu_ward.frm` | Squared Ward identity (`e1→q1` on `N` only) | prints `Sg = 0;` |
| `mumu_naive.frm` | Naive spin-summed `\|M\|^2` via `Trace4` | runs; FORM baseline |
| `mumu_schouten_c1.frm` | Appendix-C toolkit: εε norm, biline, Schouten rewrite | `NnormCheck=Biline=SchCheck=0` |
| `Makefile` | `make check` = Ward + Schouten toolkit | exit 0 |

```bash
make check          # required (Ward + Appendix-C toolkit)
make schouten       # Appendix-C toolkit only
make naive          # optional inspection of S
```

## Appendix-C toolkit notes

- FORM's Dimension-4 `e_*e_` normalization is `+24`; the Minkowski
  `(+---)` helper in Python uses an overall minus relative to FORM for
  the two-epsilon biline (`+2*(…)` in FORM vs `-2*(…)` in Python).
- The Schouten (C1) gate installs a rewrite procedure and checks it
  applies (`SchCheck=0`). It does **not** independently re-derive the
  full (C2)→(C14) chain for $M_{12}$/$M_{21}$ — that remains optional.

## What is / is not claimed here

- **Claimed:** the squared Ward premise (`Sg = 0`) on the same
  squared-trace object Section 3 studied in Python.
- **Claimed:** Appendix-C Levi-Civita normalization + biline structure +
  Schouten rewrite procedure (`make schouten`).
- **Claimed elsewhere (Python/C):** naive ≡ paper (3.5); $\sum M_{ij}/B^2$ ≡
  $L\otimes L$(naive). See report appendix “Verification of (3.5) and
  $M_{ij}$”.
- **Not claimed here:** full FORM transcription of the $M_{ij}$ blocks
  via the complete Appendix-C reduction chain.
