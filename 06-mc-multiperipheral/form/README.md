# FORM checks (muon-pair cancellation report)

Architecture choice for Section 5: **FORM**, not Python/C/Rust.
Dirac traces are what FORM is for; C owns numeric digit-loss
(Sections 6–7); the (3.5)/$M_{11}$ identities are checked in Python
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
- **Claimed elsewhere (Python):** naive ≡ paper (3.5); $M_{11}/B^2$ ≡
  $W_1\otimes W_1$(naive). See report appendix “Verification of (3.5) and
  $M_{11}$”.
- **Not claimed:** FORM transcription of $M_{12}$, $M_{21}$, $M_{22}$
  (Appendix C / beam-side $\varepsilon$ structures).
