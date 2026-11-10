# Muon-pair derivation & verification scripts

Code behind `report/mumu_cancellation.tex` Sections 3 and 5 (and the
verification appendix). Run from this directory (`scripts/`) so imports
resolve.

## Falsifiable inventory

| Script | Role | Pass criterion |
|---|---|---|
| `mumu_trace_lib.py` | Shared Dirac-trace / $N_i$ definitions | Imported; not run alone |
| `mumu_naive_trace.py` | Naive 16-term squared trace | Prints 98 monomials; writes `generated/mumu_naive_Tij.txt` |
| `mumu_gauge_check.py` | Ward identity on amplitude $N$ | Must print `ALL TERMS CANCEL` |
| `mumu_gauge_check_squared.py` | Ward identity on the squared-trace sum | Must print `VANISHES` / exit 0 |
| `mumu_stable_35_check.py` | Naive poly ≡ paper (3.5) | Must print `PASS: naive ≡ (3.5)...` |
| `mumu_mij_check.py` | (3.11) $\sum M_{ij}/B^2$ ≡ $L\otimes L$(naive), per $W_i\otimes W_j$ | Must print `PASS: M_ij/B^2...` |
| `emit_naive_poly_c.py` | Emit `numeric/eval_naive_poly.h` | Regenerates the committed C header |
| `muon_pair_denominator_check.py` | Topic 6: muon-pair vs `pi0.c` `la` bracket | Separate from this report’s Sec. 3 |

Related (not in this folder):

| Location | Role | Pass criterion |
|---|---|---|
| `form/` (`make check`) | FORM Ward + App-C + $M_{ij}$ transcription | `PASS: Sg=0` + schouten + appC + mij |
| `kinc2-driver/` (`make check`) | 2→4 gamgam/peripp + ΣMij digit-loss | `PASS: check_eemm` + `digit_loss_eemm` |
| `numeric/` (`make check`) | C (3.5) + $\sum M_{ij}$ vs Python refs | `PASS: check_stable_c` |
| topic `CMakeLists.txt` + GHA | multi-platform CMake/CTest | `ctest` 4 gates on ubuntu/macOS/Windows |
| `numeric/` (`make run`) | Sec 6 digit-loss ray naive vs (3.5) | naive digits fall; (3.5) $|S|$ keeps falling |
| `numeric/` (`make mc`) | Sec 7 MC fractions | phys `frac_bad=0`; gauge naive → 1 |

Detailed transcription notes, the $Q=p_6-p_7$ convention, the (3.5)
scanned-sign fix, the $M_{12}$/$M_{21}$ three-index $\varepsilon$
convention, and the $M_{22}$ bilinear index order are in the report
appendix “Verification of (3.5) and $M_{ij}$”.

## Run order

```bash
# Section 3 Ward + naive dump
python3 mumu_gauge_check.py
python3 mumu_naive_trace.py
python3 mumu_gauge_check_squared.py

# Section 5 follow-ups
python3 mumu_stable_35_check.py
python3 mumu_mij_check.py

# FORM + numeric + 2→4 kinc2 (from repo topic dir)
make -C ../form check
make -C ../numeric run
make -C ../numeric mc
make -C ../kinc2-driver check
```

If you change `mumu_trace_lib.py`, re-run the Section 3 trio and both
Section 5 checkers; regenerate `generated/mumu_naive_Tij.txt` before
recompiling the report appendix that includes it; run `make -C numeric regen`
only if the committed C header must track the new poly.

## Assumptions these scripts do *not* hide

- Dot products of `{e1,e2,k1,k2,p6,p7}` are independent symbols in the
  naive expansion (momentum conservation is not substituted into the
  98-term sum).
- `D1`, `D2` stay symbolic in `mumu_naive_trace.py`; only the gauge
  scripts substitute `D_i = k_i^2 - m^2`.
- On-shell `p6.p6 = p7.p7 = m^2` **is** imposed in the naive sum and in
  the squared Ward check. Dropping it makes the squared Ward check fail.
- In (3.5) / $M_{ij}$ checks, $Q = p_6 - p_7$ (not $q_1+q_2$).
- The last two bilinears of the implemented (3.5) use an overall minus
  relative to a naive scan reading (least-squares vs naive poly).
- $M_{22}$ mix bilinears are $\varepsilon(p_2,q_2;q_1,Q)$ /
  $\varepsilon(p_1,q_1;q_2,Q)$, not against $p_i$ in the second slot.
