# mc-eqc-1d

Monte-Carlo integration of a 1-D function over the interval [0,1]
using strata sampling and modifying strata sizes to make starata
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