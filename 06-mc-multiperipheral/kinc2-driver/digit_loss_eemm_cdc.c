/*
 * Soft-|t| digit loss of the CDC-60 naive matrix element on gamgam
 * samples.
 *
 *   dig_mij  — Σ M_ij/B^2, IEEE double vs long double
 *   dig_cdc  — L⊗L(naive) in cdc_fp vs the widest IEEE evaluation
 *              of the same Hessian (long double when it is wider
 *              than double, otherwise IEEE double)
 *
 * Bins match mc_digit_loss_eemm_gamgam: soft tmax<1, mid [1,10),
 * hard >=10. A point is bad when digits < 6. PASS when the soft
 * bin has CDC-naive frac_bad high and Σ M_ij frac_bad low.
 *
 * There is no long-double skip. The CDC word is 48 bits on every
 * host; the reference is IEEE, not a second CDC format. When long
 * double is not wider than double, dig_mij is reported as healthy
 * (no extended compare exists) and dig_cdc is CDC vs IEEE double.
 *
 * Build:  make -C kinc2-driver digit-loss-eemm-cdc
 * Pass:   PASS: digit_loss_eemm_cdc
 */
#include "axolib.h"
#include "eval_mij.h"
#include "eval_naive_cdc.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#ifndef N_TRIALS
#define N_TRIALS 400
#endif
#ifndef BAD_DIGITS
#define BAD_DIGITS 6.0
#endif
#ifndef N_SOFT_NEED
#define N_SOFT_NEED 80
#endif
#ifndef SOFT_FRAC_MIJ_MAX
#define SOFT_FRAC_MIJ_MAX 0.10
#endif
#ifndef SOFT_FRAC_CDC_MIN
#define SOFT_FRAC_CDC_MIN 0.35
#endif

#define TSOFT 1.0
#define TMID 10.0

enum { BIN_SOFT = 0, BIN_MID = 1, BIN_HARD = 2, N_BINS = 3 };

EEMMINPUT eemminput;

static const double me = 0.000511;
static const double mmu = 0.105658;
static const double ebeam = 15.0;

static uint64_t rng_state = 0xCDC06011ULL;

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

static void fill_rannums(double *r)
{
	int i;
	for (i = 0; i < 7; i++) {
		r[i] = urand();
		if (r[i] < 1e-12)
			r[i] = 1e-12;
		if (r[i] > 1.0 - 1e-12)
			r[i] = 1.0 - 1e-12;
	}
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

static double digits_agree(double a, double b)
{
	double rel;
	if (!isfinite(a) || !isfinite(b))
		return 0.0;
	rel = fabs(a - b) / fmax(fabs(b), 1e-300);
	return (rel > 0.0) ? -log10(rel) : 99.0;
}

static int bin_of(double tmax)
{
	if (!(tmax > 0.0) || !isfinite(tmax))
		return -1;
	if (tmax < TSOFT)
		return BIN_SOFT;
	if (tmax < TMID)
		return BIN_MID;
	return BIN_HARD;
}

int main(void)
{
	static const char *bin_name[N_BINS] = {"soft", "mid", "hard"};
	double rannums[7];
	int n_ok = 0, n_kin = 0, i, b;
	int n_bin[N_BINS], n_bad_mij[N_BINS], n_bad_cdc[N_BINS];
	double sum_dig_mij[N_BINS], sum_dig_cdc[N_BINS];
	int ld_wider = LDBL_MANT_DIG > DBL_MANT_DIG;
	double frac_mij, frac_cdc;

	memset(n_bin, 0, sizeof n_bin);
	memset(n_bad_mij, 0, sizeof n_bad_mij);
	memset(n_bad_cdc, 0, sizeof n_bad_cdc);
	memset(sum_dig_mij, 0, sizeof sum_dig_mij);
	memset(sum_dig_cdc, 0, sizeof sum_dig_cdc);

	eemminput.me = me;
	eemminput.mu = mmu;
	eemminput.ebeam = ebeam;
	eemminput.sq = 2.0 * ebeam;

	printf("# digit_loss_eemm_cdc: CDC-naive vs IEEE, ΣMij health\n");
	printf("# ebeam=%.3f  N_TRIALS=%d  bad_if_digits<%.1f  "
	       "LDBL_MANT_DIG=%d\n",
	       ebeam, N_TRIALS, BAD_DIGITS, (int)LDBL_MANT_DIG);
	printf("# PASS: soft frac_bad_cdc > %.2f  and  frac_bad_mij < %.2f"
	       "  (n_soft>=%d)\n",
	       SOFT_FRAC_CDC_MIN, SOFT_FRAC_MIJ_MAX, N_SOFT_NEED);
	printf("# bin  n  frac_bad_mij  frac_bad_cdc  mean_dig_mij  "
	       "mean_dig_cdc\n");

	for (i = 0; i < N_TRIALS; i++) {
		double dj, tmax, Sm, Sc, Sref, dig_m, dig_c;
		long double Sml;
		V4 p1, p2, p3, p5, p6, p7, q1, q2, Q;
		V4l q1l, q2l, Ql, p1l, p2l;

		fill_rannums(rannums);
		dj = gamgam(ebeam, me, me, me, me, mmu, mmu, 0.0, 2.0 * ebeam,
			    rannums, 0, 1);
		if (dj == 0.0) {
			n_kin++;
			continue;
		}
		reconstruct(&p1, &p2, &p3, &p5, &p6, &p7);
		vsub(&q1, &p1, &p3);
		vsub(&q2, &p2, &p5);
		vsub(&Q, &p6, &p7);
		tmax = fmax(fabs(extra.t1), fabs(extra.t2));
		b = bin_of(tmax);
		if (b < 0)
			continue;

		Sm = eval_mij_sum_over_B2(&q1, &q2, &Q, &p1, &p2, mmu);
		Sc = dressed_naive_cdc(&q1, &q2, &p1, &p2, &p6, &p7, mmu);
		if (ld_wider) {
			v4_to_l(&q1l, &q1);
			v4_to_l(&q2l, &q2);
			v4_to_l(&Ql, &Q);
			v4_to_l(&p1l, &p1);
			v4_to_l(&p2l, &p2);
			Sml = eval_mij_sum_over_B2_l(&q1l, &q2l, &Ql, &p1l, &p2l,
						     (long double)mmu);
			Sref = dressed_naive_ld(&q1, &q2, &p1, &p2, &p6, &p7,
						mmu);
			dig_m = digits_agree(Sm, (double)Sml);
		} else {
			Sref = dressed_naive_dbl(&q1, &q2, &p1, &p2, &p6, &p7,
						 mmu);
			dig_m = 99.0;
		}
		if (!isfinite(Sm) || !isfinite(Sc) || !isfinite(Sref))
			continue;
		dig_c = digits_agree(Sc, Sref);

		n_bin[b]++;
		sum_dig_mij[b] += dig_m;
		sum_dig_cdc[b] += dig_c;
		if (dig_m < BAD_DIGITS)
			n_bad_mij[b]++;
		if (dig_c < BAD_DIGITS)
			n_bad_cdc[b]++;
		n_ok++;
	}

	for (b = 0; b < N_BINS; b++) {
		int n = n_bin[b];
		if (n == 0) {
			printf("%s  %d  —\n", bin_name[b], n);
			continue;
		}
		printf("%s  %d  %.4f  %.4f  %.2f  %.2f\n", bin_name[b], n,
		       (double)n_bad_mij[b] / (double)n,
		       (double)n_bad_cdc[b] / (double)n,
		       sum_dig_mij[b] / (double)n, sum_dig_cdc[b] / (double)n);
	}
	printf("# accepted=%d  kin_fail=%d  n_soft=%d\n", n_ok, n_kin,
	       n_bin[BIN_SOFT]);

	if (n_ok < 100 || n_bin[BIN_SOFT] < N_SOFT_NEED) {
		fprintf(stderr, "FAIL: too few points (ok=%d soft=%d)\n", n_ok,
			n_bin[BIN_SOFT]);
		return 1;
	}
	frac_mij = (double)n_bad_mij[BIN_SOFT] / (double)n_bin[BIN_SOFT];
	frac_cdc = (double)n_bad_cdc[BIN_SOFT] / (double)n_bin[BIN_SOFT];
	if (frac_mij > SOFT_FRAC_MIJ_MAX) {
		fprintf(stderr, "FAIL: soft frac_bad_mij=%.4f > %.2f\n",
			frac_mij, SOFT_FRAC_MIJ_MAX);
		return 1;
	}
	if (frac_cdc < SOFT_FRAC_CDC_MIN) {
		fprintf(stderr, "FAIL: soft frac_bad_cdc=%.4f < %.2f\n",
			frac_cdc, SOFT_FRAC_CDC_MIN);
		return 1;
	}
	printf("PASS: digit_loss_eemm_cdc\n");
	return 0;
}
