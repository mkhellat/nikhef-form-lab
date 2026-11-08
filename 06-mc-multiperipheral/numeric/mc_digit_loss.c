/*
 * Monte Carlo digit-loss fraction for the naive 98-term poly
 * (report/mumu_cancellation.tex Section 7).
 *
 * Two scans, same random muon / polarization angles:
 *
 *   phys × eps   — transverse e1,e2; bin by photon virtuality eps
 *                  (q1^2=q2^2=-eps). Multiperipheral peak = eps→0.
 *                  Expected: frac_bad ≈ 0 (peak alone does not break
 *                  this subsystem poly).
 *
 *   gauge × delta — e1 = q1 + delta * e_T (e_T random transverse);
 *                  fixed mild eps. Bin by delta→0 (Ward surface).
 *                  Expected: frac_bad → 1 (Section 6 mechanism,
 *                  now as a phase-space fraction).
 *
 * "Bad" = digits_agree(double, long double) < BAD_DIGITS.
 *
 * Build: make mc && ./mc_digit_loss
 */

#include "eval_naive_poly.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef N_PER_BIN
#define N_PER_BIN 5000
#endif
#ifndef BAD_DIGITS
#define BAD_DIGITS 6.0
#endif

typedef struct { double e[4]; } vec4;

static uint64_t rng_state = 0xC0FFEE1234ULL;

static uint64_t xorshift64(void)
{
	uint64_t x = rng_state;
	x ^= x << 13;
	x ^= x >> 7;
	x ^= x << 17;
	rng_state = x;
	return x;
}

static double urand(void)
{
	return (xorshift64() >> 11) * (1.0 / 9007199254740992.0);
}

static double vdot(const vec4 *a, const vec4 *b)
{
	return a->e[0] * b->e[0] - a->e[1] * b->e[1] - a->e[2] * b->e[2] -
	       a->e[3] * b->e[3];
}

static void vscale(vec4 *o, const vec4 *a, double s)
{
	o->e[0] = a->e[0] * s;
	o->e[1] = a->e[1] * s;
	o->e[2] = a->e[2] * s;
	o->e[3] = a->e[3] * s;
}

static void vadd(vec4 *o, const vec4 *a, const vec4 *b)
{
	o->e[0] = a->e[0] + b->e[0];
	o->e[1] = a->e[1] + b->e[1];
	o->e[2] = a->e[2] + b->e[2];
	o->e[3] = a->e[3] + b->e[3];
}

static void fill_from_vectors(const vec4 *e1, const vec4 *e2, const vec4 *k1,
			      const vec4 *k2, const vec4 *p6, const vec4 *p7,
			      double m, Dots *d, Dotsl *dl)
{
	d->m = m;
	d->D1 = vdot(k1, k1) - m * m;
	d->D2 = vdot(k2, k2) - m * m;
	d->e1_e1 = vdot(e1, e1);
	d->e1_e2 = vdot(e1, e2);
	d->e1_k1 = vdot(e1, k1);
	d->e1_k2 = vdot(e1, k2);
	d->e1_p6 = vdot(e1, p6);
	d->e1_p7 = vdot(e1, p7);
	d->e2_e2 = vdot(e2, e2);
	d->e2_k1 = vdot(e2, k1);
	d->e2_k2 = vdot(e2, k2);
	d->e2_p6 = vdot(e2, p6);
	d->e2_p7 = vdot(e2, p7);
	d->k1_k1 = vdot(k1, k1);
	d->k1_k2 = vdot(k1, k2);
	d->k1_p6 = vdot(k1, p6);
	d->k1_p7 = vdot(k1, p7);
	d->k2_k2 = vdot(k2, k2);
	d->k2_p6 = vdot(k2, p6);
	d->k2_p7 = vdot(k2, p7);
	d->p6_p7 = vdot(p6, p7);

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

static double digits_agree(double Sd, long double Sl)
{
	double rel = fabs(Sd - (double)Sl) / fmax(fabs((double)Sl), 1e-300);
	return (rel > 0.0) ? -log10(rel) : 99.0;
}

/*
 * mode 0: physical transverse e1,e2
 * mode 1: e1 = q1 + delta * e_T  (near gauge); e2 transverse
 */
static void sample_point(double W, double eps, double delta, double m, int mode,
			 Dots *d, Dotsl *dl)
{
	double E = W / 2.0;
	double qz = sqrt(E * E + eps);
	double p = sqrt(E * E - m * m);
	double phi = 2.0 * M_PI * urand();
	double cth = 2.0 * urand() - 1.0;
	double sth = sqrt(fmax(0.0, 1.0 - cth * cth));
	double psi = 2.0 * M_PI * urand();
	double chi = 2.0 * M_PI * urand();
	vec4 q1 = {{E, 0, 0, qz}};
	vec4 q2 = {{E, 0, 0, -qz}};
	vec4 p6 = {{E, p * sth * cos(phi), p * sth * sin(phi), p * cth}};
	vec4 p7 = {{E, -p6.e[1], -p6.e[2], -p6.e[3]}};
	vec4 k1 = {{p6.e[0] - q1.e[0], p6.e[1] - q1.e[1], p6.e[2] - q1.e[2],
		    p6.e[3] - q1.e[3]}};
	vec4 k2 = {{p6.e[0] - q2.e[0], p6.e[1] - q2.e[1], p6.e[2] - q2.e[2],
		    p6.e[3] - q2.e[3]}};
	vec4 ex = {{0, 1, 0, 0}};
	vec4 ey = {{0, 0, 1, 0}};
	vec4 e1t, e2t, e1, e2, tmp;

	vscale(&e1t, &ex, cos(psi));
	vscale(&tmp, &ey, sin(psi));
	vadd(&e1t, &e1t, &tmp);

	vscale(&e2t, &ex, cos(chi));
	vscale(&tmp, &ey, sin(chi));
	vadd(&e2t, &e2t, &tmp);

	if (mode == 0) {
		e1 = e1t;
		e2 = e2t;
	} else {
		vscale(&tmp, &e1t, delta);
		vadd(&e1, &q1, &tmp);
		e2 = e2t;
	}

	fill_from_vectors(&e1, &e2, &k1, &k2, &p6, &p7, m, d, dl);
}

static void run_scan(const char *name, int mode, const double *bins, int n_bins,
		     double fixed_eps, double fixed_delta, double W, double m,
		     int n)
{
	Dots d;
	Dotsl dl;

	printf("# ensemble=%s  mode=%d  N=%d  bad_if_digits<%.1f\n", name, mode,
	       n, BAD_DIGITS);
	if (mode == 0)
		printf("# (phys) bin=eps  fixed_delta unused; e1,e2 transverse\n");
	else
		printf("# (gauge) bin=delta  fixed_eps=%.3e; e1=q1+delta*e_T\n",
		       fixed_eps);
	printf("# bin  n_bad  frac_bad  mean_digits  min_digits\n");

	for (int ib = 0; ib < n_bins; ib++) {
		double bin = bins[ib];
		double eps = (mode == 0) ? bin : fixed_eps;
		double delta = (mode == 0) ? fixed_delta : bin;
		int n_bad = 0;
		double sum_digits = 0.0;
		double min_digits = 99.0;

		rng_state = 0xC0FFEE1234ULL ^
			    (uint64_t)(ib + 1 + 17 * mode) * 0x9E3779B97F4A7C15ULL;

		for (int i = 0; i < n; i++) {
			double Sd, dig;
			long double Sl;

			sample_point(W, eps, delta, m, mode, &d, &dl);
			Sd = eval_naive_poly(&d);
			Sl = eval_naive_poly_l(&dl);
			dig = digits_agree(Sd, Sl);
			sum_digits += dig;
			if (dig < min_digits)
				min_digits = dig;
			if (dig < BAD_DIGITS)
				n_bad++;
		}
		printf("%.3e  %d  %.4f  %.2f  %.2f\n", bin, n_bad,
		       (double)n_bad / (double)n, sum_digits / (double)n,
		       min_digits);
	}
	printf("#\n");
}

int main(void)
{
	const double W = 20.0;
	const double m = 0.105658;
	const int n = N_PER_BIN;
	static const double eps_bins[] = {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6,
					  1e-7, 1e-8, 1e-9, 1e-10, 1e-11,
					  1e-12};
	static const double delta_bins[] = {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6,
					    1e-7, 1e-8, 1e-9, 1e-10};
	const int n_eps = (int)(sizeof eps_bins / sizeof eps_bins[0]);
	const int n_delta = (int)(sizeof delta_bins / sizeof delta_bins[0]);

	printf("# mc_digit_loss: fraction of naive-poly points with digit loss\n");
	printf("# W=%.3f  m=%.6f  (q1^2=q2^2=-eps)\n", W, m);

	run_scan("phys_vs_eps", 0, eps_bins, n_eps, 0.0, 0.0, W, m, n);
	run_scan("gauge_vs_delta", 1, delta_bins, n_delta, 1e-2, 0.0, W, m, n);
	return 0;
}
