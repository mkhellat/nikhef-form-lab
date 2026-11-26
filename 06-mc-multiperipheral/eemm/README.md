# eemm — VEGAS driver for e⁻e⁺ → e⁻e⁺μ⁻μ⁺

Thin replacement for the vendored `eemm.c` that `kinc2/README.md`
describes and does not ship: `vegas` on the 7-dimensional unit cube,
`gamgam` for the Jacobian, then a matrix-element backend.

Beam energy is 15 GeV by default (the same PETRA-scale value as
`kinc2-driver/check_eemm_gamgam`); an optional second argument sets
another `Ebeam`; `--ebeam GEV` or `--sqrts GEV` (= 2 `Ebeam`) do the
same. Masses are the physical electron and muon. `m4` runs
from 0 to `2 Ebeam`. There is **no** minimum-|t| cut and no
scattering-angle veto, so the soft multiperipheral peak stays inside
the integral.

| --backend | |M|² |
|---|---|
| `stable` (default) | `¼ Σ M_ij / (t1 t2 B)²`, `B = civita.bb` |
| `peripp` | `kinc2/peripp.c`, the archive's own matrix element |
| `ieee` | same `¼/(t1 t2)²` prefactor on L⊗L(naive) in IEEE double |
| `cdc` | same prefactor on L⊗L(naive) in `cdc_fp` |

The integrand is `(ħc)² * dj * |M|²` in nb, `(ħc)² = 0.389379e6 nb GeV²`.
`dj` already holds the flux `1/(2s)`, the couplings `e⁸` and the
`(2π)⁻⁸` of four-body phase space, all through `gamgam`'s `mconstant`
(the derivation is at the top of `eemm_vegas.c`). An earlier version
also divided by `2s`, copying `eee.c`; that counted the flux twice.

At `Ebeam = 15` the stable and peripp backends are gated against the
paper value 119.555 ± 0.019 nb (twophoton.pdf, PDF p.10) within
max(3 sd, 1%). Both give the same σ, about 119.67 nb: about 0.1%
high, several quoted sd, cause not pinned down. The archive run at
`Ebeam = 150` (ee.tex:63, 299822.63 ± 284 pb) is reproduced to 1e-4,
the size of our newer ħc and muon mass.

```bash
make -C 06-mc-multiperipheral/eemm check
./eemm_vegas peripp 150
./eemm_vegas ieee
./eemm_vegas cdc
./eemm_vegas --backend peripp --ebeam 150
./eemm_vegas --seed 7 --calls 50000 --iters 5
./eemm_vegas --help
```

The positional form (`BACKEND [EBEAM]`) still works and overrides
`--backend` and `--ebeam`. Other options: `--calls N`, `--iters N`,
`--acc X`. `--seed N` (1 to 2147483562; 0 is rejected) seeds the
`kinc2` `ranf` once before `vegas`, so the stream is reproducible.
Without `--seed` the fixed default stream is used and the `make check`
numbers are unchanged. The 119.555 nb gate applies only to the
reference configuration; any other run prints `# no reference gate`.
An invalid option exits with status 2.

`make check` prints `PASS: eemm_vegas` and `PASS: eemm_vegas_peripp`
(about 3.4 s each). The stable run also checks, on one accepted point
with both `|t|` in `(0.05, 5)`, that the IEEE-naive and CDC-naive
values are finite. The optional `ieee` and `cdc` runs use fewer calls
because the Hessian is much heavier than `Σ M_ij`, and have no gate.
