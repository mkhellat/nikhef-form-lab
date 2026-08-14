#include "axolib.h"

/*
	Mapping according to ds/sqrt(la(s,y,z))
	This is only useful when la(smin,y,z) can become rather small
*/

double mapla(double smin, double smax, double y, double z
					, double *ds, double x)
{
	double s,c,c1,c2,c3,xmb,xpb,alp,alm,am,ap,yy,zz,ax;
	xmb = smin-y-z;
	xpb = smax-y-z;
	c   = -4*y*z;
	c1  = xpb*xpb+c;
	c2  = xmb*xmb+c;
	if ( c1 >= 0 && c2 >= 0 ) {
		alp = sqrt(c1);
		alm = sqrt(c2);
		am  = xmb+alm;
		ap  = xpb+alp;
		yy  = ap/am;
		zz  = pow(yy,x);
		s   = y+z+(am*zz-c/(am*zz))*0.5;
		c3	= s-y-z;
		c3  = c3*c3+c;
		if ( c3 >= 0 && yy >= 0 ) {
			ax  = sqrt(c3);
			*ds = ax*log(yy);
			return(s);
		}
	}
	*ds = 0;
	return(0);
}

