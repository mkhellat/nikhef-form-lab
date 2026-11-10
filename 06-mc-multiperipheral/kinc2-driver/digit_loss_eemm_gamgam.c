/*
	Digit-loss on real 2→4 multiperipheral kinematics (gamgam):

	  Σ M_ij / B²   vs   L⊗L(naive poly)

	Samples e⁻e⁺→e⁻e⁺μ⁻μ⁺ via vendored gamgam.c, reconstructs lab
	momenta, bins by tmax = max(|t₁|,|t₂|), and for each accepted
	point reports double-vs-long-double digits of:

	  dig_mij      — eval_mij_sum_over_B2 (stable (3.11) path)
	  dig_dressed  — L⊗L contraction of eval_naive_poly (Hessian)

	plus their mutual relative agreement. Soft |t| is the
	multiperipheral peak; expect dig_mij to stay healthy while the
	dressed-naive path loses digits (same moral as Secs 6–7, now on
	real pickin/orient/gamgam phase space).

	Build:  make digit-loss-eemm
	Pass:   prints PASS: digit_loss_eemm_gamgam
*/

#include "axolib.h"
#include "eval_mij.h"
#include "eval_naive_poly.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef N_TRIALS
#define N_TRIALS 3000
#endif
#ifndef BAD_DIGITS
#define BAD_DIGITS 6.0
#endif
#ifndef N_SOFT_NEED
#define N_SOFT_NEED 80
#endif

EEMMINPUT eemminput;

static const double me = 0.000511;
static const double mmu = 0.105658;
static const double ebeam = 15.0;

static uint64_t rng_state = 0xD19170EEULL;

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

static void reconstruct_momenta(V4 *p1, V4 *p2, V4 *p3, V4 *p5, V4 *p6,
				V4 *p7)
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

static double digits_agree(double a, long double b)
{
	double rel = fabs(a - (double)b) / fmax(fabs((double)b), 1e-300);
	return (rel > 0.0) ? -log10(rel) : 99.0;
}

static void fill_dots(const V4 *e1, const V4 *e2, const V4 *q1, const V4 *q2,
		      const V4 *p6, const V4 *p7, double m, Dots *d, Dotsl *dl)
{
	V4 k1, k2;
	vsub(&k1, p6, q1);
	vsub(&k2, p6, q2);
	d->m = m;
	d->D1 = v4dot(&k1, &k1) - m * m;
	d->D2 = v4dot(&k2, &k2) - m * m;
	d->e1_e1 = v4dot(e1, e1);
	d->e1_e2 = v4dot(e1, e2);
	d->e1_k1 = v4dot(e1, &k1);
	d->e1_k2 = v4dot(e1, &k2);
	d->e1_p6 = v4dot(e1, p6);
	d->e1_p7 = v4dot(e1, p7);
	d->e2_e2 = v4dot(e2, e2);
	d->e2_k1 = v4dot(e2, &k1);
	d->e2_k2 = v4dot(e2, &k2);
	d->e2_p6 = v4dot(e2, p6);
	d->e2_p7 = v4dot(e2, p7);
	d->k1_k1 = v4dot(&k1, &k1);
	d->k1_k2 = v4dot(&k1, &k2);
	d->k1_p6 = v4dot(&k1, p6);
	d->k1_p7 = v4dot(&k1, p7);
	d->k2_k2 = v4dot(&k2, &k2);
	d->k2_p6 = v4dot(&k2, p6);
	d->k2_p7 = v4dot(&k2, p7);
	d->p6_p7 = v4dot(p6, p7);

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

typedef struct {
	const V4 *q1, *q2, *p1, *p2, *p6, *p7;
	double m;
	int use_ld;
	V4 e2_fixed;
} DressCtx;

static double S_at(const DressCtx *ctx, const V4 *e1, const V4 *e2)
{
	Dots d;
	Dotsl dl;
	fill_dots(e1, e2, ctx->q1, ctx->q2, ctx->p6, ctx->p7, ctx->m, &d, &dl);
	if (ctx->use_ld)
		return (double)eval_naive_poly_l(&dl);
	return eval_naive_poly(&d);
}

static void L_full(const V4 *q, const V4 *p, double L[4][4])
{
	double qq = v4dot(q, q);
	double pq = v4dot(p, q);
	double r[4];
	double eta[4] = {1.0, -1.0, -1.0, -1.0};
	int mu, nu;
	for (mu = 0; mu < 4; mu++)
		r[mu] = p->e[mu] - q->e[mu] * (pq / qq);
	for (mu = 0; mu < 4; mu++) {
		for (nu = 0; nu < 4; nu++) {
			double Lw1 = 2.0 * qq * ((mu == nu) ? eta[mu] : 0.0) -
				     2.0 * q->e[mu] * q->e[nu];
			double Lw2 = 8.0 * r[mu] * r[nu];
			L[mu][nu] = Lw1 + Lw2;
		}
	}
}

static double contract_L(const V4 *q, const V4 *p, DressCtx *ctx,
			 double (*f_of_e)(DressCtx *, const V4 *))
{
	double L[4][4], H[4][4];
	V4 em, en, esum;
	int mu, nu, i;
	double sum = 0.0;

	L_full(q, p, L);
	for (mu = 0; mu < 4; mu++) {
		for (i = 0; i < 4; i++)
			em.e[i] = (i == mu) ? 1.0 : 0.0;
		for (nu = mu; nu < 4; nu++) {
			for (i = 0; i < 4; i++)
				en.e[i] = (i == nu) ? 1.0 : 0.0;
			if (mu == nu) {
				H[mu][nu] = f_of_e(ctx, &em);
			} else {
				for (i = 0; i < 4; i++)
					esum.e[i] = em.e[i] + en.e[i];
				H[mu][nu] = 0.5 * (f_of_e(ctx, &esum) -
						   f_of_e(ctx, &em) -
						   f_of_e(ctx, &en));
				H[nu][mu] = H[mu][nu];
			}
		}
	}
	for (mu = 0; mu < 4; mu++)
		for (nu = 0; nu < 4; nu++)
			sum += H[mu][nu] * L[mu][nu];
	return sum;
}

static double f_e1_fixed_e2(DressCtx *ctx, const V4 *e1)
{
	return S_at(ctx, e1, &ctx->e2_fixed);
}

static double f_e2_outer(DressCtx *ctx, const V4 *e2)
{
	ctx->e2_fixed = *e2;
	return contract_L(ctx->q1, ctx->p1, ctx, f_e1_fixed_e2);
}

static double dressed_naive(const V4 *q1, const V4 *q2, const V4 *p1,
			    const V4 *p2, const V4 *p6, const V4 *p7, double m,
			    int use_ld)
{
	DressCtx ctx;
	ctx.q1 = q1;
	ctx.q2 = q2;
	ctx.p1 = p1;
	ctx.p2 = p2;
	ctx.p6 = p6;
	ctx.p7 = p7;
	ctx.m = m;
	ctx.use_ld = use_ld;
	return contract_L(q2, p2, &ctx, f_e2_outer);
}

/* log10(tmax) bins: [-12,-10), …, [-2,0) */
#define N_BINS 6

static int bin_of(double tmax)
{
	double lt;
	int b;
	if (tmax <= 0.0)
		return -1;
	lt = log10(tmax);
	if (lt < -12.0)
		lt = -12.0;
	if (lt >= 0.0)
		return N_BINS - 1;
	b = (int)floor((lt + 12.0) / 2.0);
	if (b < 0)
		b = 0;
	if (b >= N_BINS)
		b = N_BINS - 1;
	return b;
}

int main(void)
{
	double rannums[7];
	int n_ok = 0, n_kin = 0, i, b;
	int n_bin[N_BINS], n_bad_mij[N_BINS], n_bad_dr[N_BINS];
	double sum_dig_mij[N_BINS], sum_dig_dr[N_BINS], sum_agree[N_BINS];
	double soft_dig_mij = 0.0, soft_dig_dr = 0.0;
	int n_soft = 0;

	memset(n_bin, 0, sizeof n_bin);
	memset(n_bad_mij, 0, sizeof n_bad_mij);
	memset(n_bad_dr, 0, sizeof n_bad_dr);
	memset(sum_dig_mij, 0, sizeof sum_dig_mij);
	memset(sum_dig_dr, 0, sizeof sum_dig_dr);
	memset(sum_agree, 0, sizeof sum_agree);

	eemminput.me = me;
	eemminput.mu = mmu;
	eemminput.ebeam = ebeam;
	eemminput.sq = 2.0 * ebeam;

	printf("# digit_loss_eemm_gamgam: ΣMij/B^2 vs L⊗L(naive) on gamgam\n");
	printf("# ebeam=%.3f  mmu=%.6f  N_TRIALS=%d  bad_if_digits<%.1f\n",
	       ebeam, mmu, N_TRIALS, BAD_DIGITS);

	/* MSVC (and some ABIs) make long double the same width as double;
	   then dig_dressed cannot demonstrate collapse vs dig_mij. */
	if (LDBL_MANT_DIG <= DBL_MANT_DIG) {
		printf("# SKIP: long double not wider than double "
		       "(LDBL_MANT_DIG=%d DBL_MANT_DIG=%d)\n",
		       (int)LDBL_MANT_DIG, (int)DBL_MANT_DIG);
		printf("PASS: digit_loss_eemm_gamgam (skipped LD compare)\n");
		return 0;
	}
	printf("# bin_lo  bin_hi  n  frac_bad_mij  frac_bad_dr  "
	       "mean_dig_mij  mean_dig_dr  mean_agree_dig\n");

	for (i = 0; i < N_TRIALS; i++) {
		double dj, tmax, Sm, Sml, Sd, Sdl, dig_m, dig_d, dig_ag;
		V4 p1, p2, p3, p5, p6, p7, q1, q2, Q;
		V4l q1l, q2l, Ql, p1l, p2l;

		fill_rannums(rannums);
		dj = gamgam(ebeam, me, me, me, me, mmu, mmu, 0.0, 2.0 * ebeam,
			    rannums, 0, 1);
		if (dj == 0.0) {
			n_kin++;
			continue;
		}
		reconstruct_momenta(&p1, &p2, &p3, &p5, &p6, &p7);
		vsub(&q1, &p1, &p3);
		vsub(&q2, &p2, &p5);
		vsub(&Q, &p6, &p7);
		tmax = fmax(fabs(extra.t1), fabs(extra.t2));
		b = bin_of(tmax);
		if (b < 0)
			continue;

		v4_to_l(&q1l, &q1);
		v4_to_l(&q2l, &q2);
		v4_to_l(&Ql, &Q);
		v4_to_l(&p1l, &p1);
		v4_to_l(&p2l, &p2);

		Sm = eval_mij_sum_over_B2(&q1, &q2, &Q, &p1, &p2, mmu);
		Sml = eval_mij_sum_over_B2_l(&q1l, &q2l, &Ql, &p1l, &p2l,
					     (long double)mmu);
		Sd = dressed_naive(&q1, &q2, &p1, &p2, &p6, &p7, mmu, 0);
		Sdl = dressed_naive(&q1, &q2, &p1, &p2, &p6, &p7, mmu, 1);

		if (!isfinite(Sm) || !isfinite(Sd) || !isfinite((double)Sml) ||
		    !isfinite((double)Sdl))
			continue;

		dig_m = digits_agree(Sm, Sml);
		dig_d = digits_agree(Sd, (long double)Sdl);
		dig_ag = digits_agree(Sm, (long double)Sd);

		n_bin[b]++;
		sum_dig_mij[b] += dig_m;
		sum_dig_dr[b] += dig_d;
		sum_agree[b] += dig_ag;
		if (dig_m < BAD_DIGITS)
			n_bad_mij[b]++;
		if (dig_d < BAD_DIGITS)
			n_bad_dr[b]++;
		n_ok++;

		if (tmax < 1.0) {
			n_soft++;
			soft_dig_mij += dig_m;
			soft_dig_dr += dig_d;
		}
	}

	for (b = 0; b < N_BINS; b++) {
		double lo = -12.0 + 2.0 * b;
		double hi = lo + 2.0;
		int n = n_bin[b];
		if (n == 0) {
			printf("%.0e  %.0e  %d  —\n", pow(10.0, lo),
			       pow(10.0, hi), n);
			continue;
		}
		printf("%.0e  %.0e  %d  %.4f  %.4f  %.2f  %.2f  %.2f\n",
		       pow(10.0, lo), pow(10.0, hi), n,
		       (double)n_bad_mij[b] / (double)n,
		       (double)n_bad_dr[b] / (double)n,
		       sum_dig_mij[b] / (double)n, sum_dig_dr[b] / (double)n,
		       sum_agree[b] / (double)n);
	}

	printf("# accepted=%d  kin_fail=%d  soft_|t|<1: n=%d  "
	       "mean_dig_mij=%.2f  mean_dig_dressed=%.2f\n",
	       n_ok, n_kin, n_soft,
	       n_soft ? soft_dig_mij / n_soft : 0.0,
	       n_soft ? soft_dig_dr / n_soft : 0.0);

	if (n_ok < 200) {
		fprintf(stderr, "FAIL: too few accepted (%d)\n", n_ok);
		return 1;
	}
	if (n_soft < N_SOFT_NEED) {
		fprintf(stderr, "FAIL: too few soft points (%d < %d)\n", n_soft,
			N_SOFT_NEED);
		return 1;
	}
	if (soft_dig_mij / n_soft < 8.0) {
		fprintf(stderr, "FAIL: soft mean dig_mij=%.2f < 8\n",
			soft_dig_mij / n_soft);
		return 1;
	}
	/* Dressed-naive must be the worse path near the multiperipheral peak. */
	if (soft_dig_dr / n_soft >= soft_dig_mij / n_soft - 0.5) {
		fprintf(stderr,
			"FAIL: expected dig_dressed < dig_mij on soft |t| "
			"(dr=%.2f mij=%.2f)\n",
			soft_dig_dr / n_soft, soft_dig_mij / n_soft);
		return 1;
	}

	printf("PASS: digit_loss_eemm_gamgam\n");
	return 0;
}
