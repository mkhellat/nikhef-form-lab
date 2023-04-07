#include "axolib.h"

extern double s;
extern double emass;
extern double overal;
extern double Ebeam;
extern double anglecut;
extern double pi;

double meeee(double);

double feeeea(double *x)
/*
	Two dimensional kinematics function for e-e+ -> e-e+
	This time with a mapping according to
		1/(x[0]+a)^2
	This corresponds to the 1/t^2 pole which manifests itself as
	1/(1-costheta)^2 -> 1/x[0]^2.
	We cannot map that exactly so we put a small parameter a with it.
*/
{
	PVEC p_in, pa, p1, p2;
	double jacobian, t, total, st1, st2;
	double y[2],b,a,dxdu;
	p_in.e = 2*Ebeam;      pa.e = Ebeam;
	p_in.px = 0.;          pa.px = 0;
	p_in.py = 0.;          pa.py = 0;
	p_in.pz = 0.;          pa.pz = sqrt((Ebeam-emass)*(Ebeam+emass));
	p_in.m = sqrt(s);      pa.m = emass;
	p1.m = p2.m = emass;
/*
	Mapping. Our angle cut is in terms of sin(theta) but in mgoto2
	we integrate in terms of cos(theta). We translate that to a value
	for y[0] = (1-cos(theta))/2
	Then we make a mapping for the integral over y[0] from a to 1-a
	with the distribution 1/(y[0]*y[0])
*/
	a = anglecut*anglecut;
	a = a/(2*(1+sqrt(1-a)));
	b = 1/a-1/(1-a);
	y[0] = 1/(1/a-x[0]*b);
	dxdu = y[0]*y[0]*b;

	y[1] = x[1];
	if ( ( jacobian = mgoto2(&p_in,&p1,&p2,y) ) == 0 ) return(0);

	lplot.xvalue[0] = 1-2*y[0];
	lplot.xvalue[1] = 2*pi*x[1];
	lplot.xvalue[2] = x[0];
/*
	Put now cuts if needed
*/
	st1 = sqrt((p1.px*p1.px+p1.py*p1.py)/(p1.px*p1.px+p1.py*p1.py+p1.pz*p1.pz));
	st2 = sqrt((p2.px*p2.px+p2.py*p2.py)/(p2.px*p2.px+p2.py*p2.py+p2.pz*p2.pz));

	if ( ( st1 < anglecut ) || ( st2 < anglecut ) ) return(0);

	t = 2*emass*emass-2*DOTPR(pa,p1);
	total = dxdu*overal*jacobian*meeee(t);

	return(total);
}
