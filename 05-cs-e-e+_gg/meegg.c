#include "axolib.h"

extern double s;
extern double emass;

double meegg(double t, double u)
/*
		Matrix element for e-e+ -> gamma gamma
		Calculated with ex4.frm
*/
{
	double me2 = emass * emass;
	double me4 = me2*me2;
	double pe_p1 = (me2-t)/2;
	double pp_p1 = (me2-u)/2;
	double pe_p12 = pe_p1*pe_p1;
	double pp_p12 = pp_p1*pp_p1;
	return(-(24/pe_p12*me4-48/pe_p1/pp_p1*me4-8/pe_p1*pp_p1-48/
      pe_p1*me2-8*pe_p1/pp_p1+24/pp_p12*me4-48/pp_p1*me2));
}
