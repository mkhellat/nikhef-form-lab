# vegas — thin VEGAS integrator

G. P. Lepage, *J. Comput. Phys.* 27 (1978) 192, in the form shipped with
the FORM course as `formcourse/kinc/vegas.c`. The same text is already
in this lab under `02S-sample-vegas-pipeline/vegas.c` and the topic
03–05 copies. Histogram calls (`currentplot`, `Replot`, `xplot`,
`PlotIt`) are left out. Uniform random numbers and the integer/real
powers are `kinc2/ranf.c`, `ipow.c`, and `iipow.c` from the vendored
course library, not a new generator.

`avgi`, `sd`, and `chi2a` are the accumulated integral, its standard
deviation, and chi-squared per iteration, as in `kinc2/axolib.h`. The
prototype matches that header.

```bash
make -C 06-mc-multiperipheral/vegas check
```

prints `PASS: vegas_toy`. The smoke integrates `4/(1+x^2)` on
`[0,1]` (exact value π) and `exp(x+y)` on the unit square (exact
value `(e-1)^2`). A perfectly flat integrand makes an iteration
variance zero, so the course code's `1/variance` weight is
undefined; the smoke does not use one.
