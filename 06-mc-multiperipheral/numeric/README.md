# Numeric digit-loss scans (report Sections 6–7)

C evaluation of the naive 98-term squared-amplitude polynomial from
Section 3. Full PASS inventory for the report (including FORM/Python
gates) lives in `scripts/README-mumu.md` and the report verification
appendix.

| File | Role |
|---|---|
| `eval_naive_poly.h` | Auto-generated double / long-double poly (committed) |
| `digit_loss_naive.c` | Sec 6: fixed near-gauge δ scan |
| `mc_digit_loss.c` | Sec 7: MC fraction bad vs eps (phys) and vs δ (gauge) |
| `Makefile` | `make` / `make run` / `make mc` / `make regen` |

## Build and run

```bash
make
make run    # Section 6 ray
make mc     # Section 7 ensembles (N=5000 per bin)
```

### Section 6 (`make run`)

Stress: $e_1 = q_1 + \delta\,\hat x$ at fixed mild virtuality. Exact Ward
$\Rightarrow S\to0$ as $\delta\to0$. Expect digits-agree to fall from
$\sim11$ at $\delta=10^{-1}$ toward disagreement for $\delta\lesssim10^{-6}$.

### Section 7 (`make mc`)

Marks a point **bad** when `digits_agree < 6`. Expect:

- `phys_vs_eps`: `frac_bad = 0` down to `eps=1e-12`
- `gauge_vs_delta`: `frac_bad` rises toward 1 as `delta ≲ 1e-4`

## Regenerating the header

```bash
make regen   # runs scripts/emit_naive_poly_c.py
```

Do this only after changing `mumu_trace_lib.py`.
