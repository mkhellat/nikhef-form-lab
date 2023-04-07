#include "axolib.h"

extern double s;
extern double emass;

double meeee(double t)
/*
		Matrix element for e-e+ -> e-e+
		Calculated with ex2.frm
*/
{
	double s2 = s*s;
	double t2 = t*t;
	double me2 = emass * emass;
	double me4 = me2*me2;
    return(48+16/s2*t2+32/s*t+32*s/t+16*s2/t2-64*me2/s2*t
      -64*me2*s/t2+64*me4/s2-64*me4/s/t+64*me4/t2);
}
