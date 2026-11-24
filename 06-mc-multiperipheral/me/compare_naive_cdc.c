/*
 * Mild-kinematics smoke: long-double L⊗L(naive) vs the CDC-60
 * polynomial. Away from the soft-|t| peak the 48-bit coefficient
 * should still track extended precision. Pass: PASS: compare_naive_cdc
 *
 * The point is not a gamgam sample. Virtualities are O(1) GeV^2 so
 * this does not claim the soft-peak digit-loss gate (that is a later
 * driver). Build: make -C 06-mc-multiperipheral/me check
 */
#include "eval_naive_cdc.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static void setv(V4 *v, double e, double x, double y, double z)
{
	v->e[0] = e;
	v->e[1] = x;
	v->e[2] = y;
	v->e[3] = z;
}

static int close_rel(double a, double b, double rel, const char *what)
{
	double scale = fmax(fabs(a), fabs(b));
	double err;

	if (!isfinite(a) || !isfinite(b) || scale == 0.0) {
		fprintf(stderr, "FAIL: %s non-finite or zero a=%.6g b=%.6g\n",
			what, a, b);
		return 1;
	}
	err = fabs(a - b) / scale;
	printf("# %s  ld=%.8g  cdc=%.8g  rel=%.3g\n", what, b, a, err);
	if (err > rel) {
		fprintf(stderr, "FAIL: %s rel %.3g > %.3g\n", what, err, rel);
		return 1;
	}
	return 0;
}

int main(void)
{
	Dots d;
	Dotsl dl;
	V4 q1, q2, p1, p2, p6, p7;
	double poly_c, poly_l, dress_c, dress_l;
	double m = 0.105658;
	int bad = 0;

	memset(&d, 0, sizeof(d));
	d.D1 = 2.0;
	d.D2 = 3.0;
	d.m = 0.5;
	d.e1_e1 = 0.2;
	d.e1_e2 = 0.15;
	d.e1_k1 = 0.4;
	d.e1_k2 = -0.25;
	d.e1_p6 = 0.7;
	d.e1_p7 = -0.3;
	d.e2_e2 = -0.2;
	d.e2_k1 = 0.35;
	d.e2_k2 = 0.45;
	d.e2_p6 = -0.55;
	d.e2_p7 = 0.6;
	d.k1_k1 = 1.8;
	d.k1_k2 = 0.9;
	d.k1_p6 = 1.1;
	d.k1_p7 = 0.4;
	d.k2_k2 = 2.2;
	d.k2_p6 = 0.8;
	d.k2_p7 = -0.5;
	d.p6_p7 = 1.4;
	dl.D1 = d.D1;
	dl.D2 = d.D2;
	dl.m = d.m;
	dl.e1_e1 = d.e1_e1;
	dl.e1_e2 = d.e1_e2;
	dl.e1_k1 = d.e1_k1;
	dl.e1_k2 = d.e1_k2;
	dl.e1_p6 = d.e1_p6;
	dl.e1_p7 = d.e1_p7;
	dl.e2_e2 = d.e2_e2;
	dl.e2_k1 = d.e2_k1;
	dl.e2_k2 = d.e2_k2;
	dl.e2_p6 = d.e2_p6;
	dl.e2_p7 = d.e2_p7;
	dl.k1_k1 = d.k1_k1;
	dl.k1_k2 = d.k1_k2;
	dl.k1_p6 = d.k1_p6;
	dl.k1_p7 = d.k1_p7;
	dl.k2_k2 = d.k2_k2;
	dl.k2_p6 = d.k2_p6;
	dl.k2_p7 = d.k2_p7;
	dl.p6_p7 = d.p6_p7;

	poly_c = eval_naive_poly_cdc(&d);
	poly_l = (double)eval_naive_poly_l(&dl);
	/* Hand dots are O(1), not a soft peak. 1e-8 leaves margin under
	   a 48-bit coefficient (~1e-14) if the sum is only mildly ill
	   conditioned, and still fails a dropped-sign transcription. */
	bad |= close_rel(poly_c, poly_l, 1e-8, "poly");

	/* |q1^2| and |q2^2| are several GeV^2, not the |t|->0 peak. */
	setv(&q1, 4.0, 0.5, 0.2, 0.3);
	setv(&q2, 3.5, -0.4, 0.25, -0.2);
	setv(&p1, 6.0, 0.1, 0.0, 5.5);
	setv(&p2, 6.0, -0.1, 0.0, -5.5);
	setv(&p6, 2.0, 0.6, 0.3, 0.4);
	setv(&p7, 2.2, -0.3, 0.2, 0.5);
	dress_c = dressed_naive_cdc(&q1, &q2, &p1, &p2, &p6, &p7, m);
	dress_l = dressed_naive_ld(&q1, &q2, &p1, &p2, &p6, &p7, m);
	bad |= close_rel(dress_c, dress_l, 1e-8, "dressed");

	if (bad)
		return 1;
	printf("PASS: compare_naive_cdc\n");
	return 0;
}
