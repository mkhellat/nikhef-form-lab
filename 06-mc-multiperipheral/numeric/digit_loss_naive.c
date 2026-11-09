/*
 * Digit-loss scan for the naive 98-term poly vs stable Vermaseren (3.5)
 * (report/mumu_cancellation.tex Section 6).
 *
 * Stress case: photon polarization e1 = q1 + delta * n_perp (near gauge).
 * Exact Ward identity => S -> 0 as delta -> 0. The expanded naive poly
 * loses digits; the Levi-Civita form (3.5) should keep them.
 *
 * Columns:
 *   digits_naive  — agree(double, long double) of the 98-term poly
 *   digits_35     — agree(double, long double) of eval_stable_35
 *   digits_nv35   — agree(naive long double, (3.5) long double)
 *
 * Build: make -C numeric && ./digit_loss_naive
 */

#include "eval_naive_poly.h"
#include "eval_stable_35.h"

#include <math.h>
#include <stdio.h>

typedef struct { double e[4]; } vec4;

static double vdot(const vec4 *a, const vec4 *b)
{
	return a->e[0] * b->e[0] - a->e[1] * b->e[1] - a->e[2] * b->e[2] -
	       a->e[3] * b->e[3];
}

static double digits_agree_dd(double a, double b)
{
	double rel = fabs(a - b) / fmax(fabs(b), 1e-300);
	return (rel > 0.0) ? -log10(rel) : 99.0;
}

static double digits_agree_dl(double a, long double b)
{
	return digits_agree_dd(a, (double)b);
}

static void fill_all(double W, double eps, double delta, double m, Dots *d,
		     Dotsl *dl, V4 *q1, V4 *q2, V4 *e1, V4 *e2, V4 *Q,
		     V4l *q1l, V4l *q2l, V4l *e1l, V4l *e2l, V4l *Ql)
{
	double E = W / 2.0;
	double qz = sqrt(E * E + eps);
	double p = sqrt(E * E - m * m);
	vec4 vq1 = {{E, 0, 0, qz}};
	vec4 vq2 = {{E, 0, 0, -qz}};
	vec4 vp6 = {{E, p, 0, 0}};
	vec4 vp7 = {{E, -p, 0, 0}};
	vec4 vk1 = {{vp6.e[0] - vq1.e[0], vp6.e[1] - vq1.e[1],
		     vp6.e[2] - vq1.e[2], vp6.e[3] - vq1.e[3]}};
	vec4 vk2 = {{vp6.e[0] - vq2.e[0], vp6.e[1] - vq2.e[1],
		     vp6.e[2] - vq2.e[2], vp6.e[3] - vq2.e[3]}};
	/* e1 = q1 + delta * xhat  (exact gauge at delta=0) */
	vec4 ve1 = {{vq1.e[0], vq1.e[1] + delta, vq1.e[2], vq1.e[3]}};
	vec4 ve2 = {{0, 0, 1, 0}};
	vec4 vQ = {{vp6.e[0] - vp7.e[0], vp6.e[1] - vp7.e[1],
		    vp6.e[2] - vp7.e[2], vp6.e[3] - vp7.e[3]}};

	d->m = m;
	d->D1 = vdot(&vk1, &vk1) - m * m;
	d->D2 = vdot(&vk2, &vk2) - m * m;
	d->e1_e1 = vdot(&ve1, &ve1);
	d->e1_e2 = vdot(&ve1, &ve2);
	d->e1_k1 = vdot(&ve1, &vk1);
	d->e1_k2 = vdot(&ve1, &vk2);
	d->e1_p6 = vdot(&ve1, &vp6);
	d->e1_p7 = vdot(&ve1, &vp7);
	d->e2_e2 = vdot(&ve2, &ve2);
	d->e2_k1 = vdot(&ve2, &vk1);
	d->e2_k2 = vdot(&ve2, &vk2);
	d->e2_p6 = vdot(&ve2, &vp6);
	d->e2_p7 = vdot(&ve2, &vp7);
	d->k1_k1 = vdot(&vk1, &vk1);
	d->k1_k2 = vdot(&vk1, &vk2);
	d->k1_p6 = vdot(&vk1, &vp6);
	d->k1_p7 = vdot(&vk1, &vp7);
	d->k2_k2 = vdot(&vk2, &vk2);
	d->k2_p6 = vdot(&vk2, &vp6);
	d->k2_p7 = vdot(&vk2, &vp7);
	d->p6_p7 = vdot(&vp6, &vp7);

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

	for (int i = 0; i < 4; i++) {
		q1->e[i] = vq1.e[i];
		q2->e[i] = vq2.e[i];
		e1->e[i] = ve1.e[i];
		e2->e[i] = ve2.e[i];
		Q->e[i] = vQ.e[i];
	}
	v4_to_l(q1l, q1);
	v4_to_l(q2l, q2);
	v4_to_l(e1l, e1);
	v4_to_l(e2l, e2);
	v4_to_l(Ql, Q);
}

int main(void)
{
	const double W = 20.0;
	const double m = 0.105658;
	const double eps = 1e-2; /* fixed mild virtuality */
	double delta = 1e-1;
	Dots d;
	Dotsl dl;
	V4 q1, q2, e1, e2, Q;
	V4l q1l, q2l, e1l, e2l, Ql;

	printf("# digit_loss_naive: near-gauge scan naive poly vs stable (3.5)\n");
	printf("# W=%.3f  m=%.6f  eps=%.3e  (e1 = q1 + delta*xhat)\n", W, m,
	       eps);
	printf("# delta  S_naive_d  S_35_d  dig_naive  dig_35  dig_nv35\n");

	for (int i = 0; i < 12; i++) {
		double Sd, S35d, dig_n, dig_35, dig_nv;
		long double Sl, S35l;

		fill_all(W, eps, delta, m, &d, &dl, &q1, &q2, &e1, &e2, &Q,
			 &q1l, &q2l, &e1l, &e2l, &Ql);
		Sd = eval_naive_poly(&d);
		Sl = eval_naive_poly_l(&dl);
		S35d = eval_stable_35(&q1, &q2, &e1, &e2, &Q);
		S35l = eval_stable_35_l(&q1l, &q2l, &e1l, &e2l, &Ql);
		dig_n = digits_agree_dl(Sd, Sl);
		dig_35 = digits_agree_dl(S35d, S35l);
		dig_nv = digits_agree_dl((double)Sl, S35l);
		printf("%.3e  %.10e  %.10e  %.1f  %.1f  %.1f\n", delta, Sd, S35d,
		       dig_n, dig_35, dig_nv);
		delta *= 0.1;
	}
	printf("# Expect: dig_naive falls with delta; dig_35 stays high; "
	       "dig_nv35 tracks dig_naive once the poly is damaged.\n");
	return 0;
}
