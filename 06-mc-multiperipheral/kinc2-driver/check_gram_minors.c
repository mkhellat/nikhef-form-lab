/*
	Numerical confirmation (and refutation, in one specific respect) of
	Vermaseren 1983 p.353's claim, quoted in notes/topic4-pseudoscalar-
	coupling.md: "the first term [of formula (3.2)] ... is the Gram
	determinant of the system and the other 3 terms are minors of it."
	See scripts/gram_minors_pi0.py for the algebraic (symbolic) side of
	this check -- this driver checks the same claim numerically, against
	the REAL kinc2 pickin.c/orient.c/pi0.c code, at real on-shell
	kinematic points.

	Method: pickin.c already computes and stores the raw pairwise dot
	products (dotp.p12, dotp.p13, ... dotp.p2k1) it needs internally, in
	the DOTP struct. This driver reconstructs, from those same dot
	products (NOT from levi.gram/dd2/dd4's own stable factored forms),
	the naive 4x4 Gram determinant of (p1,q1,p2,q2) and its two relevant
	3x3 principal minors -- i.e. it evaluates the Levi-Civita-squared
	quantities the "slow, naive" way the paper warns is numerically
	dangerous -- and compares the result to pickin.c's own stable
	levi.gram, levi.dd2, levi.dd4.

	q1 = p1-p3, q2 = p2-p5 (pickin.c/orient.c's own definitions). With
	massless e-/e+ beams (m1=m2=m3=m5=0, used throughout this directory's
	verification tooling): p1.q1 = -dotp.p13, p2.q2 = -dotp.p25,
	q1.q2 = p1.p2 - p1.p5 - p2.p3 + p3.p5 = dotp.p12 - dotp.p15 -
	dotp.p23 + dotp.p35 (dotp.p1k2 = p1.q2, dotp.p2k1 = p2.q1 are already
	provided directly by pickin.c).

	Result (see ../notes/topic4-pseudoscalar-coupling.md for the writeup):
	At kinematic points where t1, t2 are not extremely close to their
	near-real-photon limit, the naive dot-product minors agree with
	pickin.c's stable levi.dd2 (dropping p1) and levi.dd4 (dropping p2)
	to ~0.1-1%, confirming the minor identity holds as more than an
	analogy. But at a near-real-photon point (t1 ~ 1e-11), the naive
	3x3-minor-for-part3 calculation catastrophically cancels: its three
	additive terms are individually ~1e-8 in size while their sum, and
	pickin.c's own levi.dd4, are both ~1e-9 -- a ~30x discrepancy between
	the naive sum and the stable value, with the WRONG SIGN. This is a
	direct, concrete demonstration of exactly the "very bad cancellations
	between the various terms" the uam19 lecture warns about (Topic 5's
	second open thread) -- caught here because pickin.c never performs
	this cancelling subtraction; it computes levi.dd2/dd4 via a
	numerically stable boundary-root factorization instead (pickin.c
	lines 175, 238), which this driver deliberately bypasses to see what
	the "naive" route would have given.

	Depends on the vendored kinc2 sources in ../kinc2/ (see
	../kinc2/README.md), exactly like scan_pi0_scaling.c.

	Build (from this directory, kinc2-driver/):
		cc -Wall -g -I../kinc2 -c ../kinc2/orient.c ../kinc2/pickin.c \
			../kinc2/mapt.c ../kinc2/mapw.c ../kinc2/mapla.c \
			../kinc2/ipow.c ../kinc2/iipow.c ../kinc2/ranf.c \
			../kinc2/pi0.c
		cc -Wall -g -I../kinc2 -c check_gram_minors.c
		cc -o check_gram_minors *.o -lm
		./check_gram_minors
*/

#include "axolib.h"
#include <stdio.h>

double pi0(int par);

/* Same rationale as scan_pi0_scaling.c: supply eemminput's storage
   directly rather than linking eee.o's VEGAS/plotting machinery. */
EEMMINPUT eemminput;

static double me = 0.000511;   /* electron mass, GeV */
static double mpi = 0.135;     /* pi0 mass, GeV */

/*
	Evaluates kinematics at (s, rannums), then prints:
	  - the naive 4x4 Gram determinant of (p1,q1,p2,q2) from raw dot
	    products, vs. pickin.c's stable levi.gram (part1's quantity);
	  - the naive 3x3 minor dropping p1 (q1,p2,q2), vs. levi.dd2
	    (part2's quantity);
	  - the naive 3x3 minor dropping p2 (p1,q1,q2), vs. levi.dd4
	    (part3's quantity), broken into its three additive terms so the
	    cancellation (or lack of it) is directly visible.
*/
static void check_point(double s, double r0, double r1, double r2, double r3)
{
	double rannums[4] = {r0, r1, r2, r3};
	double dj = orient(s, me, me, me, mpi, me, rannums, 0);
	double p1q1, p1q2, p1p2, q1p2, q1q2, p2q2, t1, t2;
	double a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p;
	double m11, m12, m13, m14, gram4;
	double minor_drop_p1, term1, term2, term3, minor_drop_p2;

	if (dj == 0) {
		printf("s=%-8g rannums=(%.3g,%.3g,%.3g,%.3g)  kinematics failed\n",
			s, r0, r1, r2, r3);
		return;
	}

	p1q1 = -dotp.p13;
	p2q2 = -dotp.p25;
	p1p2 = dotp.p12;
	q1p2 = dotp.p2k1;
	p1q2 = dotp.p1k2;
	q1q2 = dotp.p12 - dotp.p15 - dotp.p23 + dotp.p35;
	t1 = extra.t1;
	t2 = extra.t2;

	/* Naive 4x4 Gram determinant of (p1,q1,p2,q2) via cofactor expansion
	   along row 0. Matrix layout matches scripts/gram_minors_pi0.py. */
	a = 0;     b = p1q1;  c = p1p2;  d = p1q2;
	e = p1q1;  f = t1;    g = q1p2;  h = q1q2;
	i = p1p2;  j = q1p2;  k = 0;     l = p2q2;
	m = p1q2;  n = q1q2;  o = p2q2;  p = t2;
	m11 = f*(k*p-l*o) - g*(j*p-l*n) + h*(j*o-k*n);
	m12 = e*(k*p-l*o) - g*(i*p-l*m) + h*(i*o-k*m);
	m13 = e*(j*p-l*n) - f*(i*p-l*m) + h*(i*n-j*m);
	m14 = e*(j*o-k*n) - f*(i*o-k*m) + g*(i*n-j*m);
	gram4 = a*m11 - b*m12 + c*m13 - d*m14;

	/* 3x3 minor dropping p1 (rows/cols q1,p2,q2). */
	minor_drop_p1 = t1*(0*t2 - p2q2*p2q2)
		- q1p2*(q1p2*t2 - p2q2*q1q2)
		+ q1q2*(q1p2*p2q2 - 0*q1q2);

	/* 3x3 minor dropping p2 (rows/cols p1,q1,q2), shown term-by-term. */
	term1 = -p1q1*p1q1*t2;
	term2 = 2*p1q1*p1q2*q1q2;
	term3 = -p1q2*p1q2*t1;
	minor_drop_p2 = term1 + term2 + term3;

	printf("s=%-8g t1=%-13.5g t2=%-13.5g\n", s, t1, t2);
	printf("  part1: naive_gram4=%-14.6g  levi.gram=%-14.6g  ratio=%.6f\n",
		gram4, levi.gram, gram4 / levi.gram);
	printf("  part2: naive_minor=%-14.6g  levi.dd2 =%-14.6g  ratio=%.6f\n",
		minor_drop_p1, levi.dd2, minor_drop_p1 / levi.dd2);
	printf("  part3: term1=%-13.6g term2=%-13.6g term3=%-13.6g sum=%-13.6g\n",
		term1, term2, term3, minor_drop_p2);
	printf("         naive_minor=%-14.6g  levi.dd4 =%-14.6g  ratio=%.6f\n\n",
		minor_drop_p2, levi.dd4, minor_drop_p2 / levi.dd4);
}

int main(void)
{
	eemminput.mu = mpi;

	printf("Well-behaved points (t1, t2 not extremely close to 0):\n");
	printf("--------------------------------------------------------\n");
	check_point(20.0, 0.45, 0.55, 0.5, 0.5);
	check_point(100.0, 0.37, 0.61, 0.5, 0.5);
	check_point(400.0, 0.5, 0.5, 0.4, 0.6);

	printf("Near-real-photon point (t1 extremely small) --\n");
	printf("naive minor calculation for part3 catastrophically cancels:\n");
	printf("--------------------------------------------------------\n");
	check_point(20.0, 0.1, 0.9, 0.2, 0.8);

	return 0;
}
