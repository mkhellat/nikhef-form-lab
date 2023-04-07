## [cs-e-e+_e-e+](main4.c) - electron positron elastic scattering cross-section

Calculation of electron-positron scattering to electron-positron scattering
cross-section as observed by a detector which sweeps $cos(\theta^{CM}) > 0.05$.

A second program [main4a.c](main4a.c) and [feeeea.c](feeeea.c) that
takes advantage of _mapping_ to perform variance reduction is also
created to take care of $\frac{1}{t^2}$ integrations in the relevant
cross-section computation. This is mainly used to help VEGAS, more
quickly, focus on the most important regions, i.e. the main
contributions to this cross-section comes from $\frac{1}{t}$
components of the matrix element.

# Requirements

Generating plots requires `latex`, `dvips`, and `ps2pdf`.

# Run

Build the project using

````
make -f make4 plot
````

This would build prog4, prog4a and generate the relevant plots.
