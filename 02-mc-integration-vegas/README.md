## [mc-eqc-1d](mc-eqc-1d.c) - MC integration equalizing bin contributions

Monte-Carlo integration of a 1-D function over the interval [0,1]
using stratified sampling and modifying strata sizes to make starata
contributions to the integral close to each other

# Run

Compile the code using

````
cc mc-eqc-1d.c ranf.c -o mc-eqc-1d -lm
````

and simply run

````
./mc-str-1d <number_of_strata> <number_of_sample_points_per_stratum>
````

Upon successful execution, the program would run for 100 iterations
and report results comparing each result with the result of normal
single bon run.

````
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 iteration 99 : 0.97604238 +/- 0.00157414
   (no bins) : 0.92547647 +/- 0.09120032
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
````

## [mc-vegas78-1d](mc-vegas78-1d) - MC integration optimizations with VEGAS algorithm

Same as `mc-eqc-1d` plus a comulative result calculation using Peter
Lepage's VEGAS [1] [2]


# Run

Compile the code using

````
cc mc-vegas78-1d.c ranf.c -o mc-vegas78-1d -lm
````

and simply run

````
./mc-vegas78-1d <number_of_strata> <number_of_sample_points_per_stratum>
````

The program would run for 100 iterations with 'K' value set as 1000
and the comulative estimate of the integral value and its error would
be printed.

````
number of bins: 50 , number of points per bin: 50

.
.
.

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
iteration 100 : 0.99925153 +/- 0.00038313
cross check : 1.00067222 +/- 0.01807497
 --- total m -> 1044
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
************************************************************
*
* COMULATIVE RESULT : 0.99976868 +/- 0.00003858
*
************************************************************
````

***

[1] G. P. Lepage, “A New Algorithm for Adaptive Multidimensional Integration,” J. Comp. Phys. 27, 192–203 (1978)
[2] https://inspirehep.net/literature/1816703