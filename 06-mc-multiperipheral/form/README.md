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

Full algebraic identity between the naive sum and the paper’s
`M_{ij}` (III.11) is not claimed by these files yet: that needs a
faithful transcription of (III.11) into FORM. What is claimed is that
the gauge premise of the stable route holds in FORM on the same
squared-trace object Section 3 studied in Python.
