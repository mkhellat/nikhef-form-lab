/*
 * Digit-loss scan for Σ M_ij / B^2 (Vermaseren (3.11), electron-dressed).
 *
 * Fixed mild virtuality (eps); scan muon polar angle cosθ. Massless
 * electrons with 2 p_i·q_i = q_i^2. Reports double vs long-double digits
 * for eval_mij_sum_over_B2.
 *
 * Expect: dig_mij stays high across the scan (no Ward-surface collapse).
 * (A tiny-eps scan is a poor stress: make_massless_p becomes ill-conditioned
 * as q^2→0, which is a kinematics construction issue, not M_ij algebra.)
 *
 * Build: make digit_loss_mij && ./digit_loss_mij
 */

#include "eval_mij.h"

#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double digits_agree(double a, long double b)
{
	double rel = fabs(a - (double)b) / fmax(fabs((double)b), 1e-300);
	return (rel > 0.0) ? -log10(rel) : 99.0;
}

static int make_massless_p(const V4 *q, double pT, double phi, V4 *p)
{
	double q0 = q->e[0], qz = q->e[3];
	double t = v4dot(q, q);
	double a, b, c, disc, pz, E;
	int tries = 0;

	if (fabs(q0) < 1e-14)
		return -1;
	a = (qz / q0) * (qz / q0) - 1.0;
	b = (t * qz) / (q0 * q0);
	c = (t / (2.0 * q0)) * (t / (2.0 * q0)) - pT * pT;
	disc = b * b - 4.0 * a * c;
	while (disc < 0.0 && tries < 8) {
		pT *= 0.5;
		c = (t / (2.0 * q0)) * (t / (2.0 * q0)) - pT * pT;
		disc = b * b - 4.0 * a * c;
		tries++;
	}
	if (disc < 0.0 || fabs(a) < 1e-18)
		return -1;
	pz = (-b + sqrt(disc)) / (2.0 * a);
	E = (t / 2.0 + pz * qz) / q0;
	p->e[0] = E;
	p->e[1] = pT * cos(phi);
	p->e[2] = pT * sin(phi);
	p->e[3] = pz;
	return 0;
}

static int fill_kin(double W, double eps, double m, double cth, double phi,
		    double pT1, double phi1, double pT2, double phi2, V4 *q1,
		    V4 *q2, V4 *Q, V4 *p1, V4 *p2, V4l *q1l, V4l *q2l, V4l *Ql,
		    V4l *p1l, V4l *p2l)
{
	double E = W / 2.0;
	double qz = sqrt(E * E + eps);
	double pmu = sqrt(E * E - m * m);
	double sth = sqrt(fmax(0.0, 1.0 - cth * cth));
	V4 p6, p7;

	q1->e[0] = E;
	q1->e[1] = 0;
	q1->e[2] = 0;
	q1->e[3] = qz;
	q2->e[0] = E;
	q2->e[1] = 0;
	q2->e[2] = 0;
	q2->e[3] = -qz;
	p6.e[0] = E;
	p6.e[1] = pmu * sth * cos(phi);
	p6.e[2] = pmu * sth * sin(phi);
	p6.e[3] = pmu * cth;
	p7.e[0] = q1->e[0] + q2->e[0] - p6.e[0];
	p7.e[1] = q1->e[1] + q2->e[1] - p6.e[1];
	p7.e[2] = q1->e[2] + q2->e[2] - p6.e[2];
	p7.e[3] = q1->e[3] + q2->e[3] - p6.e[3];
	Q->e[0] = p6.e[0] - p7.e[0];
	Q->e[1] = p6.e[1] - p7.e[1];
	Q->e[2] = p6.e[2] - p7.e[2];
	Q->e[3] = p6.e[3] - p7.e[3];

	if (make_massless_p(q1, pT1, phi1, p1) != 0)
		return -1;
	if (make_massless_p(q2, pT2, phi2, p2) != 0)
		return -1;

	v4_to_l(q1l, q1);
	v4_to_l(q2l, q2);
	v4_to_l(Ql, Q);
	v4_to_l(p1l, p1);
	v4_to_l(p2l, p2);
	return 0;
}

int main(void)
{
	const double W = 20.0;
	const double m = 0.105658;
	const double eps = 1e-2;
	const double phi = 1.1;
	const double pT1 = 0.5, phi1 = 0.4;
	const double pT2 = 0.55, phi2 = 2.2;
	V4 q1, q2, Q, p1, p2;
	V4l q1l, q2l, Ql, p1l, p2l;
	double worst = 99.0;

	printf("# digit_loss_mij: ΣM_ij/B^2 double vs LD vs cosθ (eps=%.3e)\n",
	       eps);
	printf("# W=%.3f  m=%.6f  (massless p_i, 2p·q=q^2)\n", W, m);
	printf("# cth  S_mij_d  dig_mij\n");

	for (int i = 0; i < 17; i++) {
		double cth = -0.8 + 0.1 * i;
		double Sd, dig;
		long double Sl;

		if (fill_kin(W, eps, m, cth, phi, pT1, phi1, pT2, phi2, &q1,
			     &q2, &Q, &p1, &p2, &q1l, &q2l, &Ql, &p1l,
			     &p2l) != 0) {
			fprintf(stderr, "FAIL: kinematics at cth=%.2f\n", cth);
			return 1;
		}
		Sd = eval_mij_sum_over_B2(&q1, &q2, &Q, &p1, &p2, m);
		Sl = eval_mij_sum_over_B2_l(&q1l, &q2l, &Ql, &p1l, &p2l,
					    (long double)m);
		dig = digits_agree(Sd, Sl);
		if (dig < worst)
			worst = dig;
		printf("%+.2f  %.10e  %.1f\n", cth, Sd, dig);
	}
	printf("# worst dig_mij=%.1f (expect >= 8 across this mild scan)\n",
	       worst);
	if (worst < 8.0) {
		fprintf(stderr, "FAIL: worst dig_mij=%.1f < 8\n", worst);
		return 1;
	}
	printf("PASS: digit_loss_mij\n");
	return 0;
}
