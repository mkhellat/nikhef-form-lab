# kinc2-driver

Instrumented drivers against the vendored `../kinc2/` sources. These are
*our* programs; they call Vermaseren's `pickin`/`orient`/`gamgam`/`peripp`
rather than reimplementing them.

| Target | Role | Pass |
|---|---|---|
| `make check-eemm` | 2→4 eemm: gamgam kinematics + peripp vs ΣMᵢⱼ | `PASS: check_eemm_gamgam` |
| `make digit-loss-eemm` | ΣMᵢⱼ/B² vs L⊗L(naive) digit-loss on gamgam | `PASS: digit_loss_eemm_gamgam` |
| `make mc-digit-loss-eemm` | same samples as a frac_bad *rate* (soft/mid/hard |t|) | `PASS: mc_digit_loss_eemm_gamgam` |
| `make check` | all three eemm gates above | exit 0 |
| `scan_pi0_scaling` | Topic 5 s-scaling of π⁰ formula (3.2) | prints part1..sum table |
| `check_gram_minors` | Topic 4 naive Gram minors vs stable levi.* | demonstrates digit loss |
| `mc_xsec_s_scan` | crude MC σ(s) using `orient`+`pi0` | prints s-scan |

```bash
make check            # check-eemm + digit-loss-eemm + mc-digit-loss-eemm
make check-eemm       # kinematics + peripp ↔ ΣMij
make digit-loss-eemm  # dig_mij vs dig_dressed vs |t|
make mc-digit-loss-eemm  # frac_bad rate, soft tmax<1 / mid / hard

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

`mc_digit_loss_eemm_gamgam` is the matching *rate* gate (Section 7
style): same kinematics and eval paths, three bins (soft $t_{\max}<1$,
mid $[1,10)$, hard $\ge 10$), `frac_bad` when digits $<6$. PASS when
the soft bin has dressed `frac_bad` high and mij `frac_bad` low; it
self-skips when `long double` is not wider than `double`.