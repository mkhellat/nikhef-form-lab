# Muon-pair derivation scripts

Code behind `report/mumu_cancellation.tex` Section 3 and its appendices.
Run from this directory (`scripts/`) so imports resolve.

| Script | Role | Falsifiable claim |
|---|---|---|
| `mumu_trace_lib.py` | Shared Dirac-trace / $N_i$ definitions | Imported by the two trace scripts; not run alone |
| `mumu_naive_trace.py` | Naive 16-term squared trace | Prints 98 monomials; writes `generated/mumu_naive_Tij.txt` |
| `mumu_gauge_check.py` | Ward identity on amplitude $N$ | Must print `ALL TERMS CANCEL` |
| `mumu_gauge_check_squared.py` | Ward identity on the squared-trace sum | Must print `VANISHES` / exit 0 |
| `muon_pair_denominator_check.py` | Topic 6: muon-pair vs `pi0.c` `la` bracket | Separate from this report’s Sec. 3 |
| `emit_naive_poly_c.py` | Emit `numeric/eval_naive_poly.h` for Sec. 6 | Regenerates the committed C header |
| `mumu_stable_35_check.py` | Naive poly ≡ paper (3.5) Levi-Civita form | Must print `PASS: naive ≡ (3.5)...` |
| `mumu_mij_check.py` | (3.11) $M_{11}/B^2$ ≡ $W_1\otimes W_1$(naive) | Must print `PASS: M11/B^2...` |

## Run order (Section 3 checks)

```bash
python3 mumu_gauge_check.py
python3 mumu_naive_trace.py
python3 mumu_gauge_check_squared.py
```

If you change `mumu_trace_lib.py`, re-run all three. Regenerate
`generated/mumu_naive_Tij.txt` before recompiling the report appendix that
includes it.

## Assumptions these scripts do *not* hide

- Dot products of `{e1,e2,k1,k2,p6,p7}` are independent symbols in the
  naive expansion (momentum conservation is not substituted into the
  98-term sum).
- `D1`, `D2` stay symbolic in `mumu_naive_trace.py`; only the gauge
  scripts substitute `D_i = k_i^2 - m^2`.
- On-shell `p6.p6 = p7.p7 = m^2` **is** imposed in the naive sum and in
  the squared Ward check. Dropping it makes the squared Ward check fail.
