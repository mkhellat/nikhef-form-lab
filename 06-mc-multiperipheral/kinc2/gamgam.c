#include "axolib.h"

/*
	Routine uses orient/pickin to generate a 2->3 double t-channel
	configuration and then has the center line split into two.
	It first has to select a mass for the center line (w4), then
	it calls orient, next it has to carefully boost to the CM of w4
	keeping careful (means in a numerically stable way) track of the
	direction of the incoming particles q1 and q2, then do the two body
	phase space. If nm != 0 it uses a mapping that takes a combination of 
	forward and backward peaking into account (phi is considered trivial). 
	Eventually it boosts the whole thing back to the lab-frame.
	It evaluates all 4-momenta and in addition a number of Levi-Civita
	tensors that can be used in the matrix element(s).
	The return value is the Jacobian of the transformation.
	The array rannums should contain at least 7 random numbers. We assume
	that there is an axial symmetry along the beam axis.
	The outgoing particles in the center have the numbers 6 and 7 (they are
	typically the muons).
*/

VARIAD variad;
CIVITA civita;
DOTPS dotps;
EXT ext;
ANGU angu;
double qve[4];
static double pi = 3.14159265358979;
static double mconstant = 2.1868465e+10;

double gamgam(double ebeam, double m1, double m2, double m3, double m5,
			double m6, double m7, double m4min, double m4max,
			double *rannums, int options, int nm)
{
	double dj, m4, w67, wmin, wmax, s, dw, ecm6, pcm6, e3mp3, e1mp1, xx6;
	double eg, pg, pgx, pgy, pgz, pgp, pgg;
	double phicm6, cpcm6, spcm6, pcm6x, pcm6y, pcm6z, pc6z, pc6x, qcx, qcz;
	double h1, h2, hq, p6x, p6y, p6z, p7x, p7y, p7z, q0, qx, qy, qz;
	double b1, b2, b3, c1, c2, c3, r12, r13, r22, r23;
	dotps.w6 = m6*m6;
	dotps.w7 = m7*m7;
	w67 = dotps.w6-dotps.w7;
	wmin = m6+m7;
	if ( wmin < m4min ) wmin = m4min;
	wmin = wmin*wmin;
	variab.e = 2*ebeam;
	s = variab.e*variab.e;
	wmax = variab.e-m3-m5;
	if ( wmax > m4max ) wmax = m4max;
	wmax = wmax*wmax;
	pickzz.w4 = mapw(wmin,wmax,&dw,rannums[4],2);
	m4 = sqrt(pickzz.w4);
	variad.w = m4;
	dj = orient(s,m1,m2,m3,m4,m5,rannums,options);
	if ( dj == 0 || extra.t1 > 0 || extra.t2 > 0 ) return(0);

	ecm6 = (pickzz.w4+w67)/(2*m4);
	pcm6 = sqrt((ecm6-m6)*(ecm6+m6));
	dj = dj*dw*pcm6/(m4*mconstant*s);
	e3mp3 = pickzz.w3/(variab.e3+variab.p3);
	e1mp1 = pickzz.w1/(variab.e1+variab.p);

	eg = (pickzz.w4+extra.t1-extra.t2)/(2*m4);
	pg = sqrt(eg*eg-extra.t1);
	pgx = -variac.pp3*variab.cp3*variab.ct4
		  -variab.st4*(variac.de3-e1mp1+e3mp3+variab.p3*variac.al3);
	pgy = -variac.pp3*variab.sp3;
	pgz = m4*variac.de3/(variab.e4+variab.p4)
		  -variab.e4*variac.de3*variac.al4/m4
		  -variac.pp3*variab.cp3*variab.e4*variab.st4/m4
		  +variab.e4*variab.ct4/m4*(variab.p3*variac.al3+e3mp3-e1mp1);
	pgp = sqrt(pgx*pgx+pgy*pgy);
	pgg = sqrt(pgp*pgp+pgz*pgz);
	if ( ( pgg > pgp*0.9 ) && ( pgg > pg ) ) pg = pgg;
	ext.stg = pgp/pg;
	ext.cpg = pgx/pgp;
	ext.spg = pgy/pgp;
	if ( ext.stg > 1 ) return(0);
	ext.ctg = sqrt(1-ext.stg*ext.stg);
	if ( pgz < 0 ) ext.ctg = -ext.ctg;
	xx6 = rannums[5];
	dotps.q1dq2 = 0.5*(pickzz.w4-extra.t1-extra.t2);
	if ( nm != 0 ) {
		double amap, bmap, ymap, beta, ddd;
		amap = dotps.q1dq2;
		bmap = sqrt((amap*amap-extra.t1*extra.t2)*(1-4*dotps.w6/pickzz.w4));
		ymap = (amap+bmap)/(amap-bmap);
		beta = pow(ymap,2*xx6-1);
		if ( ( ymap <= 0 ) || ( beta == 1 ) ) return(0);
		xx6 = (amap/bmap*(beta-1)/(beta+1)+1)*0.5;
		if ( xx6 > 1 ) xx6 = 1;
		else if ( xx6 < 0 ) xx6 = 0;
		angu.ctcm6 = 1-2*xx6;
		ddd = (amap+bmap*angu.ctcm6)*(amap-bmap*angu.ctcm6)/amap/bmap*log(ymap);
		dj = dj*ddd*0.5;
	}
	angu.ctcm6 = 1-2*xx6;
	angu.stcm6 = 2*sqrt(xx6*(1-xx6));
	phicm6 = 2*pi*rannums[6];
	cpcm6 = cos(phicm6);
	spcm6 = sin(phicm6);
	pcm6x = pcm6*angu.stcm6*cpcm6;
	pcm6y = pcm6*angu.stcm6*spcm6;
	pcm6z = pcm6*angu.ctcm6;
	pc6z = ext.ctg*pcm6z-ext.stg*pcm6x;

	h1 = ext.stg*pcm6z+ext.ctg*pcm6x;
	pc6x = ext.cpg*h1-ext.spg*pcm6y;
	qcx = 2*pc6x;
	qcz = 2*pc6z;
	p6y = ext.cpg*pcm6y+ext.spg*h1;
	variad.e6 = (variab.e4*ecm6+variab.p4*pc6z)/m4;
	h2 = (variab.e4*pc6z+variab.p4*ecm6)/m4;
	p6x = variab.ct4*pc6x+variab.st4*h2;
	p6z = variab.ct4*h2-variab.st4*pc6x;
	qve[0] = variab.p4*qcz/m4;
	qve[2] = 2*p6y;
	hq = variab.e4*qcz/m4;
	qve[1] = variab.ct4*qcx+variab.st4*hq;
	qve[3] = variab.ct4*hq-variab.st4*qcx;
	variad.p6 = sqrt(variad.e6*variad.e6-dotps.w6);
	variad.e7 = variab.e4-variad.e6;
	variad.p7 = sqrt(variad.e7*variad.e7-dotps.w7);
	p7x = variac.pp4-p6x;
	p7y = -p6y;
	p7z = variab.p4*variab.ct4-p6z;
	variad.pp6 = sqrt(p6x*p6x+p6y*p6y);
	variad.pp7 = sqrt(p7x*p7x+p7y*p7y);
	variad.ct6 = p6z/variad.p6;
	variad.st6 = variad.pp6/variad.p6;
	variad.ct7 = p7z/variad.p7;
	variad.st7 = variad.pp7/variad.p7;
	variad.cp6 = p6x/variad.pp6;
	variad.sp6 = p6y/variad.pp6;
	variad.cp7 = p7x/variad.pp7;
	variad.sp7 = p7y/variad.pp7;
	dotps.q1dq = eg*(2*ecm6-m4)-2*pg*pcm6*angu.ctcm6;
	civita.bb = extra.t1*extra.t2+(pickzz.w4*angu.stcm6*angu.stcm6
			   +4*dotps.w6*angu.ctcm6*angu.ctcm6)*pg*pg;
	q0 = qve[0];
	qx = qve[1];
	qy = qve[2];
	qz = qve[3];
	c1 = (qx*variab.sp3-qy*variab.cp3)*variac.pp3;
	c2 = (qz*variab.e1-q0*variab.p)*variac.pp3;
	c3 = (pickzz.d1*variab.e1*variab.e1+2*pickzz.w1*variac.de3*variab.e1
		-pickzz.w1*variac.de3*variac.de3
		+variac.pp3*variac.pp3*variab.e1*variab.e1)
			/(variab.e3*variab.p+variab.p3*variab.ct3*variab.e1);
	b1 = (qx*variab.sp5-qy*variab.cp5)*variac.pp5;
	b2 = (qz*variab.e2+q0*variab.p)*variac.pp5;
	b3 = (pickzz.d2*variab.e2*variab.e2+2*pickzz.w2*variac.de5*variab.e2
		-pickzz.w2*variac.de5*variac.de5
		+variac.pp5*variac.pp5*variab.e2*variab.e2)
			/(variab.e2*variab.p5*variab.ct5-variab.e5*variab.p);
	r12 = c2*variab.sp3+qy*c3;
	r13 = -c2*variab.cp3-qx*c3;
	r22 = b2*variab.sp5+qy*b3;
	r23 = -b2*variab.cp5-qx*b3;
	civita.epsi = dotp.p12*c1*b1+r12*r22+r13*r23;
	civita.g5 = pickzz.w1*c1*c1+r12*r12+r13*r13;
	civita.g6 = pickzz.w2*b1*b1+r22*r22+r23*r23;
	civita.a5 = -(qx*variab.cp3+qy*variab.sp3)*variac.pp3*dotp.p1k2
		-(variab.e1*q0-variab.p*qz)*(variab.cp3*variab.cp5
			+variab.sp3*variab.sp5)*variac.pp3*variac.pp5
		+(variac.de5*qz+q0*(variab.p+variab.p5*variab.ct5))*c3;
	civita.a6 = -(qx*variab.cp5+qy*variab.sp5)*variac.pp5*dotp.p2k1
		-(variab.e2*q0+variab.p*qz)*(variab.cp3*variab.cp5
			+variab.sp3*variab.sp5)*variac.pp3*variac.pp5
		+(variac.de3*qz-q0*(variab.p-variab.p3*variab.ct3))*b3;
	return(dj);
}
