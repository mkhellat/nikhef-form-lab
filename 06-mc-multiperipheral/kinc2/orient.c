#include "axolib.h"

/*
	This routine converts the information of the routine pickin into
	laboratory variables, because pickin works only in terms of invariants.
	The structure variab contains the 4-vectors for use in cuts and
	histograms. The structure variac sets information for use in
	the routine gamgam and the matrix element. These variables make it
	possible to keep things numerically stable.
	Notation: e3, e4, e5 : energies
	          p3, p4, p5 : absolute value of 3-momenta
	          pp3 etc: perpendicular (=transverse) momentum.
	          ct3 etc: cos(theta3)
	          st3 etc: sim(theta3)
	          cp3 etc: cos(phi3)
	          sp3 etc: sim(phi3)
	          de3 = difference in energy: e1-e3
	          de5 = difference in energy: e2-e5
	The return value is the Jacobian of the transformation from the unit
	cube to the variables that pickin uses.
	If any error occurs, the return value is zero.
*/

VARIAB variab;
VARIAC variac;

double orient(double s, double m1, double m2, double m3, double m4, double m5
		,double *rannums, int options)
{
	double re, a1, rr, dj;
	if ( ( dj = pickin(s,m1,m2,m3,m4,m5,rannums,options) ) == 0 ) return(0);
	variab.e  = sqrt(s);
	re = 0.5/variab.e;
	variab.e1 = re*(s+pickzz.d5);
	variab.e2 = re*(s-pickzz.d5);
	variab.p  = re*pickzz.sl1;
	variac.de3= re*(extra.s2-pickzz.w3+pickzz.d5);
	variac.de5= re*(extra.s1-pickzz.w5-pickzz.d5);
	variab.e3 = variab.e1-variac.de3;
	variab.e4 = variac.de3+variac.de5;
	variab.e5 = variab.e2-variac.de5;
	if ( variab.e4 < m4 ) return(0);
	variab.p3 = sqrt(variab.e3*variab.e3-pickzz.w3);
	variab.p4 = sqrt((variab.e4-m4)*(variab.e4+m4));
	if ( variab.p4 == 0 ) return(0);
	variab.p5  = sqrt(variab.e5*variab.e5-pickzz.w5);
	variac.pp3 = sqrt(levi.dd1/s)/variab.p;
	variac.pp5 = sqrt(levi.dd3/s)/variab.p;
	variab.st3 = variac.pp3/variab.p3;
	variab.st5 = variac.pp5/variab.p5;
	if ( ( variab.st3 > 1 ) || ( variab.st5 > 1 ) ) return(0);
	variab.ct3 = sqrt(1-variab.st3*variab.st3);
	variab.ct5 = sqrt(1-variab.st5*variab.st5);
	if ( variab.e1*variab.e3 < dotp.p13 ) variab.ct3 = -variab.ct3;
	if ( variab.e2*variab.e5 > dotp.p25 ) variab.ct5 = -variab.ct5;
	variac.al3 = variab.st3*variab.st3/(1+variab.ct3);
	variac.be5 = variab.st5*variab.st5/(1-variab.ct5);
	if ( levi.dd5 < 0 ) return(0);
	variac.pp4 = sqrt(levi.dd5/s)/variab.p;
	variab.st4 = variac.pp4/variab.p4;
	if ( variab.st4 > 1 ) return(0);
	variab.ct4 = sqrt(1-variab.st4*variab.st4);
	if ( variab.e1*variab.e4 < dotp.p14 ) variab.ct4 = -variab.ct4;
	variac.al4 = 1-variab.ct4;
	variac.be4 = 1+variab.ct4;
	if ( variab.ct4 <  0 ) variac.be4 = variab.st4*variab.st4/variac.al4;
	if ( variab.ct4 >= 0 ) variac.al4 = variab.st4*variab.st4/variac.be4;
	rr  = sqrt(-levi.gram/s)/(variab.p*variac.pp4);
	variab.sp3 = rr/variac.pp3;
	variab.sp5 = -rr/variac.pp5;
	if ( ( fabs(variab.sp3) > 1 ) || ( fabs(variab.sp5) > 1 ) ) return(0);
	variab.cp3 = -sqrt(1-variab.sp3*variab.sp3);
	variab.cp5 = -sqrt(1-variab.sp5*variab.sp5);
	a1  = variac.pp3*variab.cp3-variac.pp5*variab.cp5;
	if ( fabs(variac.pp4+variac.pp3*variab.cp3+variab.cp5*variac.pp5)
			>= fabs(fabs(a1)-variac.pp4) ) {
		if ( a1 <  0 ) variab.cp5 = -variab.cp5;
		if ( a1 >= 0 ) variab.cp3 = -variab.cp3;
	}
	return(dj);
}
