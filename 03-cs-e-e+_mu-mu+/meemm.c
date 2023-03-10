#include "axolib.h"

extern double s;
extern double emass;
extern double mumass;

double meemm(double t)
/*
        Matrix element for e-e+ -> mu-mu+
        Calculated with ex1.frm
        We use special names for the powers of the variables because we
        want to avoid the use of the power functions of C, which is its
        weakest part. It also helps to avoid too many divisions.
*/
{
    double t2 = t*t;
    double me2 = emass * emass;
    double mm2 = mumass * mumass;
    double me4 = me2*me2;
    double mm4 = mm2*mm2;
    double si = 1./s;
    double si2 = si*si;
    double value = 8+16*si2*t2+16*si*t-32*mm2*si2*t+16*mm4*si2
                  -32*me2*si2*t+32*me2*mm2*si2+16*me4*si2;
    return(value);
}

