/*
 * PETRA-style VEGAS driver for e- e+ -> e- e+ mu- mu+.
 *
 * Rebuilds the thin eemm.c callback the kinc2 README describes:
 * vegas() on the 7-cube, gamgam() for the Jacobian, then a matrix
 * element. Beam energy is 15 GeV. There is no minimum-|t| cut and
 * no scattering-angle veto, so the soft multiperipheral peak stays
 * in the integral.
 *
 * Backends (argv[1], default stable):
 *   stable  1/4 * Σ M_ij / (t1 t2 B)^2   with B = civita.bb
 *   ieee    same prefactor on L⊗L(naive) in IEEE double
 *   cdc     same prefactor on L⊗L(naive) in CDC-60 soft-float
 *
 * The returned integrand is dj * |M|^2 / (2s), with s = 4 Ebeam^2.
 * The 1/(2s) factor is the flux convention used by the pi0 scan
 * (eee.c), not a new two-photon luminosity. dj already carries
 * gamgam's mconstant. This is a finite-sigma sanity check, not a
 * published PETRA cross section.
 *
 *   make -C 06-mc-multiperipheral/eemm check
 *   ./eemm_vegas [stable|ieee|cdc]
 * Pass (default): PASS: eemm_vegas
 */
#include "axolib.h"
#include "eval_mij.h"
#include "eval_naive_cdc.h"

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

EEMMINPUT eemminput;

enum { BE_STABLE = 0, BE_IEEE = 1, BE_CDC = 2 };

static const double me = 0.000511;
static const double mmu = 0.105658;
static const double ebeam = 15.0;

static int g_be = BE_STABLE;
static long n_acc = 0;

static uint64_t rng_state = 0xE2112024ULL;

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

static void set_vec(V4 *v, double e, double px, double py, double pz)
{
	v->e[0] = e;
	v->e[1] = px;
	v->e[2] = py;
	v->e[3] = pz;
}

static void vsub(V4 *o, const V4 *a, const V4 *b)
{
	o->e[0] = a->e[0] - b->e[0];
	o->e[1] = a->e[1] - b->e[1];
	o->e[2] = a->e[2] - b->e[2];
	o->e[3] = a->e[3] - b->e[3];
}

static void reconstruct(V4 *p1, V4 *p2, V4 *p3, V4 *p5, V4 *p6, V4 *p7)
{
	set_vec(p1, variab.e1, 0.0, 0.0, variab.p);
	set_vec(p2, variab.e2, 0.0, 0.0, -variab.p);
	set_vec(p3, variab.e3, variab.p3 * variab.st3 * variab.cp3,
		variab.p3 * variab.st3 * variab.sp3, variab.p3 * variab.ct3);
	set_vec(p5, variab.e5, variab.p5 * variab.st5 * variab.cp5,
		variab.p5 * variab.st5 * variab.sp5, variab.p5 * variab.ct5);
	set_vec(p6, variad.e6, variad.p6 * variad.st6 * variad.cp6,
		variad.p6 * variad.st6 * variad.sp6, variad.p6 * variad.ct6);
	set_vec(p7, variad.e7, variad.p7 * variad.st7 * variad.cp7,
		variad.p7 * variad.st7 * variad.sp7, variad.p7 * variad.ct7);
}

static double stable_me(const V4 *q1, const V4 *q2, const V4 *Q, const V4 *p1,
			const V4 *p2)
{
	double t1 = extra.t1, t2 = extra.t2, B = civita.bb;
	double sumM, den;

	if (!(t1 < 0.0) || !(t2 < 0.0) || !isfinite(B) || B == 0.0)
		return NAN;
	sumM = eval_M11(q1, q2, Q, mmu) + eval_M12(q1, q2, Q, p2, mmu) +
	       eval_M21(q1, q2, Q, p1, mmu) + eval_M22(q1, q2, Q, p1, p2);
	den = t1 * t2 * B;
	return 0.25 * sumM / (den * den);
}

static double naive_me(double dressed)
{
	double t1 = extra.t1, t2 = extra.t2, den;

	if (!isfinite(dressed) || !(t1 < 0.0) || !(t2 < 0.0))
		return NAN;
	den = t1 * t2;
	return 0.25 * dressed / (den * den);
}

static double me_of(int be, const V4 *q1, const V4 *q2, const V4 *Q,
		    const V4 *p1, const V4 *p2, const V4 *p6, const V4 *p7)
{
	if (be == BE_IEEE)
		return naive_me(dressed_naive_dbl(q1, q2, p1, p2, p6, p7, mmu));
	if (be == BE_CDC)
		return naive_me(dressed_naive_cdc(q1, q2, p1, p2, p6, p7, mmu));
	return stable_me(q1, q2, Q, p1, p2);
}

static double feemm(double *x)
{
	double dj, me2, s;
	V4 p1, p2, p3, p5, p6, p7, q1, q2, Q;

	dj = gamgam(ebeam, me, me, me, me, mmu, mmu, 0.0, 2.0 * ebeam, x, 0, 1);
	if (!(dj > 0.0) || !isfinite(dj))
		return 0.0;
	reconstruct(&p1, &p2, &p3, &p5, &p6, &p7);
	vsub(&q1, &p1, &p3);
	vsub(&q2, &p2, &p5);
	vsub(&Q, &p6, &p7);
	me2 = me_of(g_be, &q1, &q2, &Q, &p1, &p2, &p6, &p7);
	if (!isfinite(me2))
		return 0.0;
	n_acc++;
	s = 4.0 * ebeam * ebeam;
	return dj * me2 / (2.0 * s);
}

static int moderate_point(double *st, double *ie, double *cd)
{
	int n;

	for (n = 0; n < 4000; n++) {
		double r[7], dj;
		V4 p1, p2, p3, p5, p6, p7, q1, q2, Q;
		int k;

		for (k = 0; k < 7; k++) {
			r[k] = urand();
			if (r[k] < 1e-12)
				r[k] = 1e-12;
			if (r[k] > 1.0 - 1e-12)
				r[k] = 1.0 - 1e-12;
		}
		dj = gamgam(ebeam, me, me, me, me, mmu, mmu, 0.0, 2.0 * ebeam, r,
			    0, 1);
		if (!(dj > 0.0))
			continue;
		if (fabs(extra.t1) < 0.05 || fabs(extra.t1) > 5.0)
			continue;
		if (fabs(extra.t2) < 0.05 || fabs(extra.t2) > 5.0)
			continue;
		reconstruct(&p1, &p2, &p3, &p5, &p6, &p7);
		vsub(&q1, &p1, &p3);
		vsub(&q2, &p2, &p5);
		vsub(&Q, &p6, &p7);
		*st = me_of(BE_STABLE, &q1, &q2, &Q, &p1, &p2, &p6, &p7);
		*ie = me_of(BE_IEEE, &q1, &q2, &Q, &p1, &p2, &p6, &p7);
		*cd = me_of(BE_CDC, &q1, &q2, &Q, &p1, &p2, &p6, &p7);
		printf("# point |t1|=%.4g |t2|=%.4g stable=%.6g ieee=%.6g "
		       "cdc=%.6g\n",
		       fabs(extra.t1), fabs(extra.t2), *st, *ie, *cd);
		return 0;
	}
	fprintf(stderr, "FAIL: no moderate-|t| point\n");
	return 1;
}

static const char *be_name(int be)
{
	if (be == BE_IEEE)
		return "ieee";
	if (be == BE_CDC)
		return "cdc";
	return "stable";
}

int main(int argc, char **argv)
{
	long ncalls = 2000;
	int itmax = 4;
	double st, ie, cd;

	if (argc == 2) {
		if (strcmp(argv[1], "stable") == 0)
			g_be = BE_STABLE;
		else if (strcmp(argv[1], "ieee") == 0)
			g_be = BE_IEEE;
		else if (strcmp(argv[1], "cdc") == 0)
			g_be = BE_CDC;
		else {
			fprintf(stderr, "usage: %s [stable|ieee|cdc]\n", argv[0]);
			return 2;
		}
	} else if (argc != 1) {
		fprintf(stderr, "usage: %s [stable|ieee|cdc]\n", argv[0]);
		return 2;
	}

	/* CDC and IEEE Hessians are much heavier than Σ M_ij. */
	if (g_be != BE_STABLE) {
		ncalls = 400;
		itmax = 2;
	}

	eemminput.me = me;
	eemminput.mu = mmu;
	eemminput.ebeam = ebeam;
	eemminput.sq = 2.0 * ebeam;

	printf("# eemm_vegas: Ebeam=%.3f GeV backend=%s flux=1/(2s) "
	       "no min-|t| cut\n",
	       ebeam, be_name(g_be));
	vegas(feemm, 0.05, 7, ncalls, itmax, 0);
	printf("# sigma=%.8g  sd=%.3g  chi2=%.3g  n_acc=%ld\n", avgi, sd, chi2a,
	       n_acc);

	if (!isfinite(avgi) || !(avgi > 0.0) || !isfinite(sd) || n_acc < 50) {
		fprintf(stderr, "FAIL: sigma not finite and positive\n");
		return 1;
	}
	if (g_be != BE_STABLE) {
		printf("PASS: eemm_vegas_%s\n", be_name(g_be));
		return 0;
	}
	if (moderate_point(&st, &ie, &cd))
		return 1;
	if (!isfinite(st) || !(st > 0.0) || !isfinite(ie) || !isfinite(cd)) {
		fprintf(stderr, "FAIL: backend point not finite\n");
		return 1;
	}
	printf("PASS: eemm_vegas\n");
	return 0;
}
