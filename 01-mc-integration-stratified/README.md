# mc-str-1d

Monte-Carlo integration of a 1-D function over the interval [0,1]
using strata sampling.

# Run

Compile the code using

````
cc mc-str-1d.c ranf.c -o mc-str-1d -lm
````

and simply run

````
./mc-str-1d <number_of_strata> <number_of_sample_points_per_stratum>
````

Upon successful execution, a sample result would look like this

````
::::::::::::::::::::::: FINAL RESULT :::::::::::::::::::::::
: For one range of 20 points: 0.36176291 +/- 0.03454181     
: For 10 ranges of each 2 points: 0.37014151 +/- 0.00187797 
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
````
