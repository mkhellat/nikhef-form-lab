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
| `Makefile` | `make check` gates on the Ward identity | exit 0 |

```bash
make check          # required
make naive          # optional inspection of S
```

## What is / is not claimed here

- **Claimed:** the squared Ward premise (`Sg = 0`) on the same
  squared-trace object Section 3 studied in Python.
- **Claimed elsewhere (Python):** naive ≡ paper (3.5); $\sum M_{ij}/B^2$ ≡
  $L\otimes L$(naive) per $W_i\otimes W_j$. See report appendix
  “Verification of (3.5) and $M_{ij}$”.
- **Not claimed here:** FORM transcription of the $M_{ij}$ blocks
  (Appendix C / beam-side $\varepsilon$ structures) — that gate lives in
  Python.
