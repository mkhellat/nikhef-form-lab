# eemm — VEGAS driver for e⁻e⁺ → e⁻e⁺μ⁻μ⁺

Thin replacement for the vendored `eemm.c` that `kinc2/README.md`
describes and does not ship: `vegas` on the 7-dimensional unit cube,
`gamgam` for the Jacobian, then a matrix-element backend.

Beam energy is 15 GeV (the same PETRA-scale value as
`kinc2-driver/check_eemm_gamgam`). Masses are the physical electron
and muon. `m4` runs from 0 to `2 Ebeam`. There is **no** minimum-|t|
cut and no scattering-angle veto, so the soft multiperipheral peak
stays inside the integral.

| argv | |M|² |
|---|---|
| `stable` (default) | `¼ Σ M_ij / (t1 t2 B)²`, `B = civita.bb` |
| `ieee` | same `¼/(t1 t2)²` prefactor on L⊗L(naive) in IEEE double |
| `cdc` | same prefactor on L⊗L(naive) in `cdc_fp` |

The integrand is `dj * |M|² / (2s)` with `s = 4 Ebeam²`. The `1/(2s)`
factor is the flux convention of the π⁰ scan (`eee.c`), and `dj`
already includes `gamgam`'s `mconstant`. The number is a finite-σ
sanity check, not a published PETRA cross section. On soft `|t|` the
stable formula is the one `check_eemm_gamgam` compares to `peripp`.

```bash
make -C 06-mc-multiperipheral/eemm check
./eemm_vegas ieee
./eemm_vegas cdc
```

`make check` prints `PASS: eemm_vegas`. It integrates the stable
backend and, on one accepted point with both `|t|` in `(0.05, 5)`,
checks that the IEEE-naive and CDC-naive values are finite. The
optional `ieee` and `cdc` runs use fewer calls because the Hessian
is much heavier than `Σ M_ij`.
