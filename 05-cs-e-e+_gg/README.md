## [cs-e-e+_gg](main5.c) - electron positron to gamma gamma scattering cross-section

Calculation of electron-positron scattering cross-section to a pair of
photons as observed by a detector which sweeps $cos(\theta^{CM}) > 0.05$.

The corresponding matrix element expression includes $t$ and $u$,
hence [kinematic](feegg.c) should be revised accordingly as the photon
masss is zero and the defintion of $u$ has to be included.

# Requirements

Generating plots requires `latex`, `dvips`, and `ps2pdf`.

# Run

Build the project using

````
make -f make5 plot
````

This would build prog5 and generate the relevant plots (a replica of
instance 3 has been declared for instance 5).
