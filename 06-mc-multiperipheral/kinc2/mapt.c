#include "axolib.h"

double mapt(double tmin, double tmax, double *dt, double x, int type)
{
	double y, t;
	switch ( type ) {
		case 0:
/*
			Flat (cq no) mapping
*/
			*dt = tmin - tmax;
			return(tmin - x * (*dt));
		case 1:
/*
			Assumes dt/t.   t < 0
*/
			y = tmax/tmin;
			t = tmin * pow(y,x);
			*dt = -t*log(y);
			return(t);
		default:
			printf("case %d not implemented in mapt\n",type);
			exit(-1);
			return(0);
	}
}

