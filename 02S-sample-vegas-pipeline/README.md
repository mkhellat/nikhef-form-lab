## [mc-vegas-pipeline](main2.c) - sample vegas pipeline

Monte-Carlo integration of the 2-D function fun2() over the interval
[0,1] using vegas. This is a copy of the code used in the
[uam2019](https://www.nikhef.nl/~form/maindir/courses/uam2019/uam2019.html)
course in the final segment of part3 with few comments.

# Requirements

Generating plots requires `latex`, `dvips`, and `ps2pdf`.

# Run

Build the project using

````
make -f make2 plot
````

This would build prog2 and generate the plot with

````
./prog2 plot2.tex
latex plot2
dvips plot2 -o
ps2pdf plot2.ps
````

To clean the build do

````
make -f make2 clean
````
