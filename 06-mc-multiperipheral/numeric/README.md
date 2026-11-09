# Numeric digit-loss scans (report Sections 6–7)

C evaluation of the naive 98-term squared-amplitude polynomial from
Section 3, side-by-side with the stable Vermaseren (3.5) form, plus a
smoke gate for the (3.11) $\sum M_{ij}/B^2$ evaluator. Full PASS inventory
for the report (including FORM/Python gates) lives in
`scripts/README-mumu.md` and the report verification appendix.

| File | Role |
|---|---|
| `eval_naive_poly.h` | Auto-generated double / long-double poly (committed) |
| `eval_stable_35.h` | Hand C twin of (3.5); matches `mumu_stable_35_check.py` |
| `eval_mij.h` | Hand C twin of $\sum M_{ij}/B^2$; matches `mumu_mij_check.py` |
| `digit_loss_naive.c` | Sec 6: near-gauge δ scan naive vs (3.5) |
| `mc_digit_loss.c` | Sec 7: MC frac_bad for naive and (3.5) |
| `check_stable_c.c` | Fixed-kinematics gate: (3.5) and $M_{ij}$ vs Python refs |
| `Makefile` | `make` / `make run` / `make mc` / `make check` / `make regen` |

## Build and run

```bash
make
make check  # must print PASS: check_stable_c
make run    # Section 6 ray
make mc     # Section 7 ensembles (N=5000 per bin)
```

### Section 6 (`make run`)

Stress: $e_1 = q_1 + \delta\,\hat x$ at fixed mild virtuality. Exact Ward
$\Rightarrow S\to0$ as $\delta\to0$. Columns:

- `dig_naive` — double vs long double of the 98-term poly (falls)
- `dig_35` — double vs long double of (3.5)
- `dig_nv35` — naive long double vs (3.5) long double

Expect the naive poly to saturate near $\sim10^{-13}$ while (3.5) continues
toward zero.

### Section 7 (`make mc`)

Marks a point **bad** when `digits_agree < 6`. Reports `frac_bad_naive`
and `frac_bad_35`. Expect:

- `phys_vs_eps`: both fractions near 0 down to `eps=1e-12`
- `gauge_vs_delta`: `frac_bad_naive` rises toward 1 as `delta ≲ 1e-4`;
  (3.5) is healthier at the smallest $\delta$ (absolute value keeps falling)

### $M_{ij}$ in C

`eval_mij.h` evaluates $\sum M_{ij}/B^2$ (electron-dressed (3.11)). The
Sec 6/7 drivers compare the *subreaction* object $S(e_1,e_2)$ for which
(3.5) is the stable twin; $M_{ij}$ needs beam momenta $p_1,p_2$ and is
gated by `make check` against the Python reference.

## Regenerating the naive header

```bash
make regen   # runs scripts/emit_naive_poly_c.py
```

Do this only after changing `mumu_trace_lib.py`.
