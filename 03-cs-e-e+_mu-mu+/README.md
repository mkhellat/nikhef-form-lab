## [cs-e-e+_mu-mu+](main3.c) - electron pair to muon pair scattering cross-section

Calculation of electron-positron scattering to muon-antimuon
cross-section as observed by a detector which sweeps $-0.8 <
cos(\theta^{CM}) < 0.8$ and has a single dead zone $10 < \theta^{CM} <
25, 17 < \phi^{CM} < 39$.

As explained in part4 notes of
[uam2019](https://www.nikhef.nl/~form/maindir/courses/uam2019/uam2019.html),
kinematics is aminly controlled from [feemm.c](feemm.c).

# Requirements

Generating plots requires `latex`, `dvips`, and `ps2pdf`.

# Run

Build the project using

````
make -f make3 plot
````

This would build prog3 and generate the plot like in 02S.
