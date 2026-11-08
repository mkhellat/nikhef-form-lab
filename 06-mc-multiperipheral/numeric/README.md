# Numeric digit-loss scan (report Section 6)

C evaluation of the naive 98-term squared-amplitude polynomial from
Section 3, under a near-gauge stress: $e_1 = q_1 + \delta\,\hat x$.
Exact Ward identity $\Rightarrow S\to 0$ as $\delta\to 0$, while individual
monomials stay large — the cancellation Vermaseren warned about.

| File | Role |
|---|---|
| `eval_naive_poly.h` | Auto-generated double / long-double poly (committed) |
| `digit_loss_naive.c` | Near-gauge $\delta$ scan; prints digits-agree |
| `Makefile` | `make` / `make run` / `make regen` |

## Build and run

```bash
make
make run
```

Expect digits-agree to fall from $\sim 11$ at $\delta=10^{-1}$ toward
disagreement for $\delta\lesssim 10^{-6}$.

## Regenerating the header

Requires sympy and the Section 3 scripts:

```bash
make regen   # runs scripts/emit_naive_poly_c.py
```

Do this only after changing `mumu_trace_lib.py`.
