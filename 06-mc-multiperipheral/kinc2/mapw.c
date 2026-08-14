#include "axolib.h"

double mapw(double wmin, double wmax, double *dw, double x, int type)
{
	double y, w, w1, w2;
	switch ( type ) {
		case 0:
/*
			Flat (cq no) mapping
*/
			*dw = wmax - wmin;
			return(wmin + x * (*dw));
		case 1:
/*
			Assumes dw2/w2.   w2 > 0
*/
			y = wmax/wmin;
			w = wmin * pow(y,x);
			*dw = w*log(y);
			return(w);
		case 2:
/*
			Assumes dw2/(w2*w2).
*/
			w1 = 1/wmin;
			w2 = 1/wmax;
			*dw = w1-w2;
			w = 1/(w2+(*dw)*x);
			*dw *= w*w;
			return(w);
		default:
			printf("case %d not implemented in mapw\n",type);
			exit(-1);
			return(0);
	}
}

