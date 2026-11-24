/*
 * L⊗L Hessian of the naive polynomial, same contraction as
 * kinc2-driver/digit_loss_eemm_gamgam.c (S_at → contract_L on q1,p1
 * then on q2,p2). mode selects the polynomial arithmetic.
 */
#include "eval_naive_cdc.h"

#include <string.h>

enum { MODE_LD = 0, MODE_CDC = 1, MODE_DBL = 2 };

typedef struct {
	const V4 *q1, *q2, *p1, *p2, *p6, *p7;
	double m;
	int mode;
	V4 e2_fixed;
} DressCtx;

static void vsub(V4 *o, const V4 *a, const V4 *b)
{
	o->e[0] = a->e[0] - b->e[0];
	o->e[1] = a->e[1] - b->e[1];
	o->e[2] = a->e[2] - b->e[2];
	o->e[3] = a->e[3] - b->e[3];
}

static void fill_dots(const V4 *e1, const V4 *e2, const DressCtx *ctx, Dots *d,
		      Dotsl *dl)
{
	V4 k1, k2;

	vsub(&k1, ctx->p6, ctx->q1);
	vsub(&k2, ctx->p6, ctx->q2);
	d->m = ctx->m;
	d->D1 = v4dot(&k1, &k1) - ctx->m * ctx->m;
	d->D2 = v4dot(&k2, &k2) - ctx->m * ctx->m;
	d->e1_e1 = v4dot(e1, e1);
	d->e1_e2 = v4dot(e1, e2);
	d->e1_k1 = v4dot(e1, &k1);
	d->e1_k2 = v4dot(e1, &k2);
	d->e1_p6 = v4dot(e1, ctx->p6);
	d->e1_p7 = v4dot(e1, ctx->p7);
	d->e2_e2 = v4dot(e2, e2);
	d->e2_k1 = v4dot(e2, &k1);
	d->e2_k2 = v4dot(e2, &k2);
	d->e2_p6 = v4dot(e2, ctx->p6);
	d->e2_p7 = v4dot(e2, ctx->p7);
	d->k1_k1 = v4dot(&k1, &k1);
	d->k1_k2 = v4dot(&k1, &k2);
	d->k1_p6 = v4dot(&k1, ctx->p6);
	d->k1_p7 = v4dot(&k1, ctx->p7);
	d->k2_k2 = v4dot(&k2, &k2);
	d->k2_p6 = v4dot(&k2, ctx->p6);
	d->k2_p7 = v4dot(&k2, ctx->p7);
	d->p6_p7 = v4dot(ctx->p6, ctx->p7);

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

static double S_at(DressCtx *ctx, const V4 *e1, const V4 *e2)
{
	Dots d;
	Dotsl dl;

	fill_dots(e1, e2, ctx, &d, &dl);
	if (ctx->mode == MODE_CDC)
		return eval_naive_poly_cdc(&d);
	if (ctx->mode == MODE_DBL)
		return eval_naive_poly(&d);
	return (double)eval_naive_poly_l(&dl);
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

static double f_e1_fixed_e2(DressCtx *ctx, const V4 *e1)
{
	return S_at(ctx, e1, &ctx->e2_fixed);
}

static double contract_L(const V4 *q, const V4 *p, DressCtx *ctx,
			 double (*f_of_e)(DressCtx *, const V4 *));

static double f_e2_outer(DressCtx *ctx, const V4 *e2)
{
	ctx->e2_fixed = *e2;
	return contract_L(ctx->q1, ctx->p1, ctx, f_e1_fixed_e2);
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

static double dressed_mode(const V4 *q1, const V4 *q2, const V4 *p1,
			   const V4 *p2, const V4 *p6, const V4 *p7, double m,
			   int mode)
{
	DressCtx ctx;

	memset(&ctx, 0, sizeof(ctx));
	ctx.q1 = q1;
	ctx.q2 = q2;
	ctx.p1 = p1;
	ctx.p2 = p2;
	ctx.p6 = p6;
	ctx.p7 = p7;
	ctx.m = m;
	ctx.mode = mode;
	return contract_L(q2, p2, &ctx, f_e2_outer);
}

double dressed_naive_cdc(const V4 *q1, const V4 *q2, const V4 *p1, const V4 *p2,
			 const V4 *p6, const V4 *p7, double m)
{
	return dressed_mode(q1, q2, p1, p2, p6, p7, m, MODE_CDC);
}

double dressed_naive_ld(const V4 *q1, const V4 *q2, const V4 *p1, const V4 *p2,
			const V4 *p6, const V4 *p7, double m)
{
	return dressed_mode(q1, q2, p1, p2, p6, p7, m, MODE_LD);
}

double dressed_naive_dbl(const V4 *q1, const V4 *q2, const V4 *p1, const V4 *p2,
			 const V4 *p6, const V4 *p7, double m)
{
	return dressed_mode(q1, q2, p1, p2, p6, p7, m, MODE_DBL);
}
