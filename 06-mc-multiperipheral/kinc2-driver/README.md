# kinc2-driver

Instrumented drivers against the vendored `../kinc2/` sources. These are
*our* programs; they call Vermaseren's `pickin`/`orient`/`gamgam`/`peripp`
rather than reimplementing them.

| Target | Role | Pass |
|---|---|---|
| `make check-eemm` | 2→4 eemm: gamgam kinematics + peripp vs ΣMᵢⱼ | `PASS: check_eemm_gamgam` |
| `scan_pi0_scaling` | Topic 5 s-scaling of π⁰ formula (3.2) | prints part1..sum table |
| `check_gram_minors` | Topic 4 naive Gram minors vs stable levi.* | demonstrates digit loss |
| `mc_xsec_s_scan` | crude MC σ(s) using `orient`+`pi0` | prints s-scan |

```bash
make check-eemm   # required 2→4 gate
```

`check_eemm_gamgam` reconstructs lab 4-momenta from `VARIAB`/`VARIAD`,
checks conservation / on-shell / `tᵢ`, evaluates `peripp(0,0)`, and on
soft multiperipheral samples (`|t₁|,|t₂|<1`) compares that to
`¼ ΣMᵢⱼ/(t₁ t₂ B)²` from `../numeric/eval_mij.h`.
