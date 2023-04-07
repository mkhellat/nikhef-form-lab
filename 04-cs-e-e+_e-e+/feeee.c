#include "axolib.h"

extern double s;
extern double emass;
extern double overal;
extern double Ebeam;
extern double anglecut;
extern double pi;

double meeee(double);

double feeee(double *x)
/*
	Two dimensional kinematics function for e-e+ -> e-e+
*/
{
	PVEC p_in, pa, p1, p2;
	double jacobian, t, total, st1, st2;
	p_in.e = 2*Ebeam;      pa.e = Ebeam;
	p_in.px = 0.;          pa.px = 0;
	p_in.py = 0.;          pa.py = 0;
	p_in.pz = 0.;          pa.pz = sqrt((Ebeam-emass)*(Ebeam+emass));
	p_in.m = sqrt(s);      pa.m = emass;
	p1.m = p2.m = emass ;
	if ( ( jacobian = mgoto2(&p_in,&p1,&p2,x) ) == 0 ) return(0);
	t = 2*emass*emass-2*DOTPR(pa,p1);
/*
	Put now cuts if needed
*/
	st1 = sqrt((p1.px*p1.px+p1.py*p1.py)/(p1.px*p1.px+p1.py*p1.py+p1.pz*p1.pz));
	st2 = sqrt((p2.px*p2.px+p2.py*p2.py)/(p2.px*p2.px+p2.py*p2.py+p2.pz*p2.pz));
/*
	Next we fill the histograms.
	We have defined 2 histograms in PlotBoundaries.
	The boundaries are given in 'instance' = 4.
*/
	lplot.xvalue[0] = 1-2*x[0];
	lplot.xvalue[1] = 2*pi*x[1];
	lplot.xvalue[2] = x[0];

	if ( ( st1 < anglecut ) || ( st2 < anglecut ) ) return(0);
	total = overal*jacobian*meeee(t);

	return(total);
}
