# me — naive |M|² backends

`eval_naive_cdc.c` is the 98-monomial polynomial from
`numeric/eval_naive_poly.h` (itself emitted by
`scripts/emit_naive_poly_c.py`) with every product, quotient, sum, and
integer power evaluated by `cdc_fp`. `dressed_naive_cdc` is the same
L⊗L Hessian as `kinc2-driver/digit_loss_eemm_gamgam.c`: lepton tensor
and Hessian combination stay in host double; the polynomial does not.

`dressed_naive_ld` is the long-double reference through that same
contraction. `dressed_naive_dbl` is the IEEE-double polynomial
through it, for the eemm driver.

```bash
make -C 06-mc-multiperipheral/me check
```

prints `PASS: compare_naive_cdc`. The smoke uses O(1) hand dots and a
mild virtuality (|q²| of several GeV²), where the 48-bit coefficient
should track long double. It does not sample the soft-|t| peak.
