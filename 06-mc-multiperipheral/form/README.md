# FORM checks (muon-pair cancellation report)

Architecture choice for Section 5: **FORM**, not Python/C/Rust.
Dirac traces and (later) Schouten reductions are what FORM is for; C
is reserved for numeric digit-loss drivers (Section 6); Rust would not
buy domain leverage here.

| File | Role | Pass criterion |
|---|---|---|
| `mumu_ward.frm` | Squared Ward identity (`e1→q1` on `N` only) | prints `Sg = 0;` |
| `mumu_naive.frm` | Naive spin-summed `\|M\|^2` via `Trace4` | runs; used as FORM baseline |
| `Makefile` | `make check` gates on the Ward identity | exit 0 |

```bash
make check          # required
make naive          # optional inspection of S
```

The naive ≡ paper (3.5) subreaction identity is checked numerically in
`scripts/mumu_stable_35_check.py` (Levi-Civita bilinears; see report
Sec. 5). A FORM transcription of the electron-dressed `M_{ij}` (III.11)
is still open. What these FORM files claim is the squared Ward premise
(`Sg = 0`) on the same object Section 3 studied in Python.
