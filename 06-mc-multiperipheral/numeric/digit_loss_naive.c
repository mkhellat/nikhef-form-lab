/*
 * Digit-loss scan for the naive 98-term 4-vector-product polynomial
 * (report/mumu_cancellation.tex Section 6).
 *
 * Stress case: photon polarization e1 = q1 + delta * n_perp (near gauge).
 * Exact Ward identity => S -> 0 as delta -> 0, while individual monomials
 * in the expanded poly stay large — the cancellation Vermaseren warned
 * about. Compare eval_naive_poly (double) vs eval_naive_poly_l (long
 * double) at identical inputs.
 *
 * Build: make -C numeric && ./digit_loss_naive
 */

#include "eval_naive_poly.h"

#include <math.h>
#include <stdio.h>

typedef struct { double e[4]; } vec4;

static double vdot(const vec4 *a, const vec4 *b)
{
	return a->e[0] * b->e[0] - a->e[1] * b->e[1] - a->e[2] * b->e[2] -
	       a->e[3] * b->e[3];
}

static void fill_dots(double W, double eps, double delta, double m, Dots *d,
		      Dotsl *dl)
{
	double E = W / 2.0;
	double qz = sqrt(E * E + eps);
	double p = sqrt(E * E - m * m);
	vec4 q1 = {{E, 0, 0, qz}};
	vec4 q2 = {{E, 0, 0, -qz}};
	vec4 p6 = {{E, p, 0, 0}};
	vec4 p7 = {{E, -p, 0, 0}};
	vec4 k1 = {{p6.e[0] - q1.e[0], p6.e[1] - q1.e[1], p6.e[2] - q1.e[2],
		    p6.e[3] - q1.e[3]}};
	vec4 k2 = {{p6.e[0] - q2.e[0], p6.e[1] - q2.e[1], p6.e[2] - q2.e[2],
		    p6.e[3] - q2.e[3]}};
	/* e1 = q1 + delta * xhat  (exact gauge at delta=0) */
	vec4 e1 = {{q1.e[0], q1.e[1] + delta, q1.e[2], q1.e[3]}};
	vec4 e2 = {{0, 0, 1, 0}};

	d->m = m;
	d->D1 = vdot(&k1, &k1) - m * m;
	d->D2 = vdot(&k2, &k2) - m * m;
	d->e1_e1 = vdot(&e1, &e1);
	d->e1_e2 = vdot(&e1, &e2);
	d->e1_k1 = vdot(&e1, &k1);
	d->e1_k2 = vdot(&e1, &k2);
	d->e1_p6 = vdot(&e1, &p6);
	d->e1_p7 = vdot(&e1, &p7);
	d->e2_e2 = vdot(&e2, &e2);
	d->e2_k1 = vdot(&e2, &k1);
	d->e2_k2 = vdot(&e2, &k2);
	d->e2_p6 = vdot(&e2, &p6);
	d->e2_p7 = vdot(&e2, &p7);
	d->k1_k1 = vdot(&k1, &k1);
	d->k1_k2 = vdot(&k1, &k2);
	d->k1_p6 = vdot(&k1, &p6);
	d->k1_p7 = vdot(&k1, &p7);
	d->k2_k2 = vdot(&k2, &k2);
	d->k2_p6 = vdot(&k2, &p6);
	d->k2_p7 = vdot(&k2, &p7);
	d->p6_p7 = vdot(&p6, &p7);

	dl->m = (long double)d->m;
	dl->D1 = (long double)d->D1;
	dl->D2 = (long double)d->D2;
	dl->e1_e1 = (long double)d->e1_e1;
	dl->e1_e2 = (long double)d->e1_e2;
	dl->e1_k1 = (long double)d->e1_k1;
	dl->e1_k2 = (long double)d->e1_k2;
	dl->e1_p6 = (long double)d->e1_p6;
	dl->e1_p7 = (long double)d->e1_p7;
	dl->e2_e2 = (long double)d->e2_e2;
	dl->e2_k1 = (long double)d->e2_k1;
	dl->e2_k2 = (long double)d->e2_k2;
	dl->e2_p6 = (long double)d->e2_p6;
	dl->e2_p7 = (long double)d->e2_p7;
	dl->k1_k1 = (long double)d->k1_k1;
	dl->k1_k2 = (long double)d->k1_k2;
	dl->k1_p6 = (long double)d->k1_p6;
	dl->k1_p7 = (long double)d->k1_p7;
	dl->k2_k2 = (long double)d->k2_k2;
	dl->k2_p6 = (long double)d->k2_p6;
	dl->k2_p7 = (long double)d->k2_p7;
	dl->p6_p7 = (long double)d->p6_p7;
}

int main(void)
{
	const double W = 20.0;
	const double m = 0.105658;
	const double eps = 1e-2; /* fixed mild virtuality */
	double delta = 1e-1;
	Dots d;
	Dotsl dl;

	printf("# digit_loss_naive: near-gauge scan of the 98-term poly\n");
	printf("# W=%.3f  m=%.6f  eps=%.3e  (e1 = q1 + delta*xhat)\n", W, m,
	       eps);
	printf("# delta  S_double  S_long  rel_err  digits_agree\n");

	for (int i = 0; i < 12; i++) {
		double Sd, rel, digits;
		long double Sl;

		fill_dots(W, eps, delta, m, &d, &dl);
		Sd = eval_naive_poly(&d);
		Sl = eval_naive_poly_l(&dl);
		rel = fabs(Sd - (double)Sl) /
		      fmax(fabs((double)Sl), 1e-300);
		digits = (rel > 0.0) ? -log10(rel) : 99.0;
		printf("%.3e  %.10e  %.10Le  %.3e  %.1f\n", delta, Sd, Sl, rel,
		       digits);
		delta *= 0.1;
	}
	printf("# Exact gauge (delta=0) => S=0 analytically; smaller digits "
	       "agree => cancellation damage in double.\n");
	return 0;
}
