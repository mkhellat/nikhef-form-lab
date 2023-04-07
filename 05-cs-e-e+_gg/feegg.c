#include "axolib.h"

extern double s;
extern double emass;
extern double overal;
extern double Ebeam;
extern double anglecut;
extern double pi;

double meegg(double,double);

double feegg(double *x)
/*
	Two dimensional kinematics function for e-e+ -> gamma gamma
*/
{
	PVEC p_in, pa, p1, p2;
	double jacobian, t, u, total, st1, st2;
	p_in.e = 2*Ebeam;      pa.e = Ebeam;
	p_in.px = 0.;          pa.px = 0;
	p_in.py = 0.;          pa.py = 0;
	p_in.pz = 0.;          pa.pz = sqrt((Ebeam-emass)*(Ebeam+emass));
	p_in.m = sqrt(s);      pa.m = emass;
/*
        Next we need to provide the masses of the outgoing particles.
        They are needed by mgoto2.
*/
	p1.m = p2.m = 0 ;
/*
        Now the call to mgoto2. In principle it should return a nonzero
        value unless we are exactly at the edge of phase space.
*/
	if ( ( jacobian = mgoto2(&p_in,&p1,&p2,x) ) == 0 ) return(0);
/*
        Convert the values of the 4-momenta into variables that we like more
*/
	t = emass*emass-2*DOTPR(pa,p1);
	u = emass*emass-2*DOTPR(pa,p2);
/*
	Put now cuts if needed
*/
	st1 = sqrt((p1.px*p1.px+p1.py*p1.py)/(p1.px*p1.px+p1.py*p1.py+p1.pz*p1.pz));
	st2 = sqrt((p2.px*p2.px+p2.py*p2.py)/(p2.px*p2.px+p2.py*p2.py+p2.pz*p2.pz));
/*
	Next we fill the histograms.
	We have defined 2 histograms in PlotBoundaries.
	The boundaries are given in 'instance' = 5.
*/
	lplot.xvalue[0] = 1-2*x[0];
	lplot.xvalue[1] = 2*pi*x[1];
	lplot.xvalue[2] = x[0];

	if ( ( st1 < anglecut ) || ( st2 < anglecut ) ) return(0);
/*
        Evaluate the matrix element and compose the value of the output
*/
	total = overal*jacobian*meegg(t,u);
	return(total);
}
