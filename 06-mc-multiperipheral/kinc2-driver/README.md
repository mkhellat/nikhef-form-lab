# kinc2-driver

Instrumented drivers against the vendored `../kinc2/` sources. These are
*our* programs; they call Vermaseren's `pickin`/`orient`/`gamgam`/`peripp`
rather than reimplementing them.

| Target | Role | Pass |
|---|---|---|
| `make check-eemm` | 2→4 eemm: gamgam kinematics + peripp vs ΣMᵢⱼ | `PASS: check_eemm_gamgam` |
| `make digit-loss-eemm` | ΣMᵢⱼ/B² vs L⊗L(naive) digit-loss on gamgam | `PASS: digit_loss_eemm_gamgam` |
| `make check` | both eemm gates above | exit 0 |
| `scan_pi0_scaling` | Topic 5 s-scaling of π⁰ formula (3.2) | prints part1..sum table |
| `check_gram_minors` | Topic 4 naive Gram minors vs stable levi.* | demonstrates digit loss |
| `mc_xsec_s_scan` | crude MC σ(s) using `orient`+`pi0` | prints s-scan |

```bash
make check            # check-eemm + digit-loss-eemm
make check-eemm       # kinematics + peripp ↔ ΣMij
make digit-loss-eemm  # dig_mij vs dig_dressed vs |t|

# Portable / CI (from topic root):
cmake -S .. -B ../build && cmake --build ../build
ctest --test-dir ../build --output-on-failure
```

`check_eemm_gamgam` reconstructs lab 4-momenta from `VARIAB`/`VARIAD`,
checks conservation / on-shell / `tᵢ`, evaluates `peripp(0,0)`, and on
soft multiperipheral samples (`|t₁|,|t₂|<1`) compares that to
`¼ ΣMᵢⱼ/(t₁ t₂ B)²` from `../numeric/eval_mij.h`.

`digit_loss_eemm_gamgam` bins the same gamgam samples by
$\max(|t_1|,|t_2|)$ and compares double-vs-long-double digits of
$\sum M_{ij}/B^2$ against the L⊗L Hessian of the naive poly. Soft $|t|$
keeps `dig_mij` healthy while `dig_dressed` collapses — Secs 6–7's moral
on real $2\to4$ phase space.