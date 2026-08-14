#include "axolib.h"

/*
	Routine does a 2->3 phase space conversion according to the configuration

         m1  p1 ------>------------->----- p3  m3
                             |
                             v q1
                             |
                             |------>----- p4  m4
                             |
                             ^ q2
                             |
         m2  p2 ------>------------->----- p5  m5

	The variables used are in principle t1=(p1-p3)^2, t2=(p2-p5)^2,
	s2=(p4+p5)^2 and s1=(p3+p4)^2. This transformation has the Jacobian
	1/sqrt(-Delta4(p1,p2,p3,p4)) in which Delta4 is the Gram determinant
	of the system. As it turns out, numerical stability can be improved
	by using Delta4 itself as an integration variable instead of s1.
	(the computation of its boundaries requires much care).

	The order of integration between s2,t1,t2 can be determined with the
	option parameter (order is the order in which we determine its value):
		0  : s2,t1,t2
		<0 : t1,s2,t2
		>0 : t1,t2,s2
	The absolute value of the option parameter determines the type of
	mapping used for s2. If it is >1 a special mapping according to a
	lambda function is used. When it is 1, a more regular ds2/s2 mapping
	is selected.
	The array rannums should contain at least 4 numbers between 0 and 1.
	The return value is the Jacobian of the transformation.

	The routine also computes the various dotproducts and a number of
	contractions between Levi-Civita tensors. Here those can be computed
	in a numerically stable way. They allow a numerically stable evaluation
	of the matrix element(s).
*/

PICKZZ pickzz;
EXTRA extra;
ACCURA accura;
LEVI levi;
DOTP dotp;

static double pi = 3.14159265358979;

double pickin(double s,
			  double m1, double m2, double m3, double m4, double m5,
			  double *rannums, int option)
{
	double sig, sig1, sig2, d6, ss, rl1, smax, ds2, sp, d3, rl2, sl2;
	double t1max, t1min, dt1, t2max, t2min, dt2, d4, d8, t13, sl3, s2max;
	double splus, sb, sd, se, s2x, r1,r2, rl4, sl4, r3, r4, b, c, t25, sl6;
	double sl7, sl5, s2p, s2min, ap, yy4, dd, st, delb, dj, s1p;
	double sbb, sdd, see, s1m, ssb, ssd, sse, s1pm, s1pp;
	pickzz.w1 = m1*m1;
	pickzz.w2 = m2*m2;
	pickzz.w3 = m3*m3;
	pickzz.w4 = m4*m4;
	pickzz.w5 = m5*m5;
	sig = m4+m5;
	sig1 = sig*sig;
	sig2 = sig1;
	pickzz.d1 = pickzz.w3-pickzz.w1;
	pickzz.d2 = pickzz.w5-pickzz.w2;
	pickzz.d5 = pickzz.w1-pickzz.w2;
	d6 = pickzz.w4-pickzz.w5;
	ss = s+pickzz.d5;
	rl1 = ss*ss-4*pickzz.w1*s;

	if ( rl1 <= 0 ) return(0);

	pickzz.sl1 = sqrt(rl1);

	if ( option == 0 ) {
		smax = s+pickzz.w3-2*m3*sqrt(s);
		extra.s2 = mapw(sig1,smax,&ds2,rannums[2],1);
		if ( ds2 == 0 ) return(0);
		sig1 = extra.s2;
	}
	sp = s+pickzz.w3-sig1;
	d3 = sig1-pickzz.w2;
	rl2 = sp*sp-4*s*pickzz.w3;

	if ( rl2 <= 0 ) return(0);

	sl2 = sqrt(rl2);
	t1max = pickzz.w1+pickzz.w3-(ss*sp+pickzz.sl1*sl2)/(2*s);
	t1min = (pickzz.d1*d3+(d3-pickzz.d1)
				*(d3*pickzz.w1-pickzz.d1*pickzz.w2)/s)/t1max;
	extra.t1 = mapt(t1min,t1max,&dt1,rannums[0],1);
	if ( dt1 == 0 ) return(0);

	d4 = pickzz.w4-extra.t1;
	d8 = extra.t1-pickzz.w2;
	t13 = extra.t1-pickzz.w1-pickzz.w3;
	levi.sa1 = -(extra.t1-pickzz.d1)*(extra.t1-pickzz.d1)*0.25
	           +pickzz.w1*extra.t1;
	if ( levi.sa1 >= 0 ) return(0);
	sl3 = sqrt(-levi.sa1);

	if ( pickzz.w1 == 0 ) {
		s2max = (s*(extra.t1*(s+d8-pickzz.w3)-pickzz.w2*pickzz.w3)
				+pickzz.w2*pickzz.w3*(pickzz.w2+pickzz.w3-extra.t1))/ss/t13;
		splus = sig2;
	}
	else {
		sb = (s*(extra.t1-pickzz.d1)+pickzz.d5*t13)/(2*pickzz.w1)+pickzz.w3;
		sd = pickzz.sl1*sl3/pickzz.w1;
		se = (s*(extra.t1*(s+t13-pickzz.w2)-pickzz.w2*pickzz.d1)
				+pickzz.w3*(pickzz.d5*d8+pickzz.w2*pickzz.w3))/pickzz.w1;
		if ( fabs((sb-sd)/sd) < 1 ) {
			s2max = sb+sd;
			splus = se/s2max;
		}
		else {
			splus = sb-sd;
			s2max = se/splus;
		}
	}
	s2x = s2max;
	if ( option < 0 ) {
		if ( splus > sig2 ) sig2 = splus;
		if ( option < -1 ) extra.s2 = mapw(sig2,s2max,&ds2,rannums[2],1);
		if ( option == -1 ) extra.s2 = mapla(sig2,s2max,extra.t1,pickzz.w2,&ds2,rannums[2]);
		if ( ds2 == 0 ) return(0);
	}
	if ( option <= 0 ) {
		s2x = extra.s2;
	}
	r1 = s2x-d8;
	r2 = s2x-d6;
	rl4 = (r1*r1-4*pickzz.w2*s2x)*(r2*r2-4*pickzz.w5*s2x);

	if ( rl4 <= 0 ) return(0);
	sl4 = sqrt(rl4);
	t2max = pickzz.w2+pickzz.w5-(r1*r2+sl4)/(2*s2x);
	t2min = (pickzz.d2*d4+(d4-pickzz.d2)*
			(d4*pickzz.w2-pickzz.d2*extra.t1)/s2x)/t2max;
	extra.t2 = mapt(t2min,t2max,&dt2,rannums[1],1);
	if ( dt2 == 0 ) return(0);

	pickzz.d7 = extra.t1-extra.t2;
	r3 = d4-extra.t2;
	r4 = pickzz.d2-extra.t2;
	b = r3*r4-2*(extra.t1+pickzz.w2)*extra.t2;
	c = extra.t2*d6*d8+(d6-d8)*(d6*pickzz.w2-d8*pickzz.w5);
	t25 = extra.t2-pickzz.w2-pickzz.w5;
	levi.sa2 = -r4*r4*0.25+pickzz.w2*extra.t2;
	if ( levi.sa2 >= 0 ) return(0);
	sl6 = 2*sqrt(-levi.sa2);

	levi.g4 = -0.25*r3*r3+extra.t1*extra.t2;
	if ( levi.g4 >= 0 ) return(0);
	sl7 = 2*sqrt(-levi.g4);
	sl5 = sl6*sl7;
	if ( fabs((sl5-b)/sl5) > 0 ) {
		s2p = (sl5-b)/(2*extra.t2);
		s2min = c/(extra.t2*s2p);
	}
	else {
		s2min = (-sl5-b)/(2*extra.t2);
		s2p = c/(extra.t2*s2min);
	}
	if ( option > 1 ) extra.s2 = mapw(s2min,s2max,&ds2,rannums[2],1);
	if ( option == 1 ) extra.s2 = mapla(s2min,s2max,extra.t1,pickzz.w2,&ds2,rannums[2]);
	if ( ds2 == 0 ) return(0);
	ap = -(extra.s2+d8)*(extra.s2+d8)*0.25+extra.s2*extra.t1;
	if ( pickzz.w1 == 0 )
		levi.dd1 = ss*t13*(extra.s2-s2max)*0.25;
	else
		levi.dd1 = -pickzz.w1*(extra.s2-s2max)*(extra.s2-splus)*0.25;
	levi.dd2 = -extra.t2*(extra.s2-s2p)*(extra.s2-s2min)*0.25;

	yy4 = cos(pi*rannums[3]);
	dd  = levi.dd1*levi.dd2;
	dotp.p12 = 0.5*(s-pickzz.w1-pickzz.w2);
	st  = extra.s2-extra.t1-pickzz.w2;
	delb = (2*pickzz.w2*r3+r4*st)*(4*dotp.p12*extra.t1
			-(extra.t1-pickzz.d1)*st)/(16*ap);
	if ( dd <= 0 ) return(0);
	levi.delta = delb-yy4*st*sqrt(dd)/(2*ap);
	extra.s1 = extra.t2+pickzz.w1+(2*dotp.p12*r3-4*levi.delta)/st;
	if ( ap >= 0 ) return(0);
	dj   = ds2*dt1*dt2*pi*pi/(8*pickzz.sl1*sqrt(-ap));
	levi.gram = (1-yy4)*(1+yy4)*dd/ap;

	dotp.p13  = -t13*0.5;
	dotp.p14  = (pickzz.d7+extra.s1-pickzz.w3)*0.5;
	dotp.p15  = (s+extra.t2-extra.s1-pickzz.w2)*0.5;
	dotp.p23  = (s+extra.t1-extra.s2-pickzz.w1)*0.5;
	dotp.p24  = (extra.s2-pickzz.d7-pickzz.w5)*0.5;
	dotp.p25  = -t25*0.5;
	dotp.p34  = (extra.s1-pickzz.w3-pickzz.w4)*0.5;
	dotp.p35  = (s+pickzz.w4-extra.s1-extra.s2)*0.5;
	dotp.p45  = (extra.s2-pickzz.w4-pickzz.w5)*0.5;
	dotp.p1k2 = (extra.s1-extra.t2-pickzz.w1)*0.5;
	dotp.p2k1 = st*0.5;

	if ( pickzz.w2 == 0 ) {
		s1p = (s*(extra.t2*(s-pickzz.w5+extra.t2-pickzz.w1)
				-pickzz.w1*pickzz.w5)+pickzz.w1*pickzz.w5*(pickzz.w1
				+pickzz.w5-extra.t2))/t25/(s-pickzz.d5);
		levi.dd3 = -t25*(s-pickzz.d5)*(s1p-extra.s1)*0.25;
	}
	else {
		sbb = (s*(extra.t2-pickzz.d2)-pickzz.d5*t25)/(2*pickzz.w2)+pickzz.w5;
		sdd = pickzz.sl1*sl6/(2*pickzz.w2);
		see = (s*(extra.t2*(s+t25-pickzz.w1)-pickzz.w1*pickzz.d2)
				+pickzz.w5*(pickzz.w1*pickzz.w5
				-pickzz.d5*(extra.t2-pickzz.w1)))/pickzz.w2;
		if ( sbb/sdd < 0 ) {
			s1m = sbb-sdd;
			s1p = see/s1m;
		}
		else {
			s1p = sbb+sdd;
			s1m = see/s1p;
		}
		levi.dd3 = -pickzz.w2*(s1p-extra.s1)*(s1m-extra.s1)*0.25;
	}
	accura.acc3 = (s1p-extra.s1)/(s1p+extra.s1);
	ssb = extra.t2+pickzz.w1-r3*(pickzz.d1-extra.t1)/(2*extra.t1);
	ssd = sl3*sl7/extra.t1;
	sse = (extra.t2-pickzz.w1)*(pickzz.w4-pickzz.w3)
		+(extra.t2-pickzz.w4+pickzz.d1)*((extra.t2-pickzz.w1)*
			pickzz.w3-(pickzz.w4-pickzz.w3)*pickzz.w1)/extra.t1;
	if ( ssb/ssd < 0 ) {
		s1pm = ssb-ssd;
		s1pp = sse/s1pm;
	}
	else {
		s1pp = ssb+ssd;
		s1pm = sse/s1pp;
	}
	levi.dd4 = -extra.t1*(extra.s1-s1pp)*(extra.s1-s1pm)*0.25;
	accura.acc4 = (extra.s1-s1pm)/(extra.s1+s1pm);
	levi.dd5 = levi.dd1+levi.dd3+((dotp.p12*(extra.t1-pickzz.d1)*0.5
				-pickzz.w1*dotp.p2k1)*(dotp.p2k1*(extra.t2-pickzz.d2)
						-pickzz.w2*r3)
				-levi.delta*(2*dotp.p12*dotp.p2k1
						-pickzz.w2*(extra.t1-pickzz.d1)))/dotp.p2k1;
	return(dj);
}
