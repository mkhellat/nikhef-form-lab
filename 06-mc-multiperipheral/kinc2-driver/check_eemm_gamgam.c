/*
	Genuine 2→4 multiperipheral kinematics gate for
	e⁻e⁺ → e⁻e⁺μ⁻μ⁺, mirroring what check_gram_minors.c did for the
	π⁰ (2→3) path.

	Instruments the REAL kinc2 gamgam.c / orient.c / pickin.c /
	peripp.c code (vendored from kinc2.tar.gz). Does not reimplement
	the phase-space map. For each accepted point:

	  1. Reconstruct lab 4-momenta from VARIAB/VARIAD angles
	     (φ₄ = 0 frame used inside gamgam/orient).
	  2. Check on-shell masses, 4-momentum conservation, and
	     tᵢ = qᵢ² vs extra.t1/t2 (absolute+relative tolerance —
	     near-real-photon points cancel in p1−p3).
	  3. Evaluate peripp(0,0) and require finite, non-negative.
	  4. On soft multiperipheral samples (|t₁|,|t₂| < 1 GeV²),
	     compare peripp to ¼ ΣMᵢⱼ/(t₁ t₂ B)² from
	     numeric/eval_mij.h. Hard-|t| points are *not* used for
	     that compare: angle reconstruction then loses digits
	     relative to peripp's Levi-Civita path (the paper's point).

	Build (from this directory, kinc2-driver/):
		make check-eemm
*/

#include "axolib.h"
#include "eval_mij.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef N_TRIALS
#define N_TRIALS 4000
#endif
#ifndef N_SOFT_NEED
#define N_SOFT_NEED 40
#endif

EEMMINPUT eemminput;

static const double me = 0.000511;
static const double mmu = 0.105658;
static const double ebeam = 15.0;

static uint64_t rng_state = 0xE2E42026ULL;

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

static void vadd4(V4 *o, const V4 *a, const V4 *b)
{
	o->e[0] = a->e[0] + b->e[0];
	o->e[1] = a->e[1] + b->e[1];
	o->e[2] = a->e[2] + b->e[2];
	o->e[3] = a->e[3] + b->e[3];
}

static void reconstruct_momenta(V4 *p1, V4 *p2, V4 *p3, V4 *p4, V4 *p5,
				V4 *p6, V4 *p7)
{
	set_vec(p1, variab.e1, 0.0, 0.0, variab.p);
	set_vec(p2, variab.e2, 0.0, 0.0, -variab.p);
	set_vec(p3, variab.e3, variab.p3 * variab.st3 * variab.cp3,
		variab.p3 * variab.st3 * variab.sp3, variab.p3 * variab.ct3);
	set_vec(p4, variab.e4, variab.p4 * variab.st4, 0.0,
		variab.p4 * variab.ct4);
	set_vec(p5, variab.e5, variab.p5 * variab.st5 * variab.cp5,
		variab.p5 * variab.st5 * variab.sp5, variab.p5 * variab.ct5);
	set_vec(p6, variad.e6, variad.p6 * variad.st6 * variad.cp6,
		variad.p6 * variad.st6 * variad.sp6, variad.p6 * variad.ct6);
	set_vec(p7, variad.e7, variad.p7 * variad.st7 * variad.cp7,
		variad.p7 * variad.st7 * variad.sp7, variad.p7 * variad.ct7);
}

static int t_match(double t_code, double t_vec)
{
	double abserr = fabs(t_vec - t_code);
	double scale = fmax(fabs(t_code), 1e-30);
	return abserr <= 1e-8 + 1e-3 * scale;
}

static int soft_t(double t1, double t2)
{
	return fabs(t1) < 1.0 && fabs(t2) < 1.0;
}

int main(void)
{
	double rannums[7];
	int n_ok = 0, n_kin_fail = 0, n_mass_fail = 0, n_cons_fail = 0;
	int n_t_fail = 0, n_peripp_bad = 0;
	int n_soft = 0, n_soft_fail = 0;
	int i;
	double max_cons = 0.0, max_soft_rel = 0.0;

	eemminput.me = me;
	eemminput.mu = mmu;
	eemminput.ebeam = ebeam;
	eemminput.sq = 2.0 * ebeam;

	printf("check_eemm_gamgam: ebeam=%.3f GeV, N_TRIALS=%d\n", ebeam,
	       N_TRIALS);

	for (i = 0; i < N_TRIALS; i++) {
		double dj, per, t1, t2, Bcode, sumM, pred, rel;
		double m2_3, m2_5, m2_6, m2_7, cons;
		double t1v, t2v;
		V4 p1, p2, p3, p4, p5, p6, p7, q1, q2, Q, sum_in, sum_out, tmp;

		fill_rannums(rannums);
		dj = gamgam(ebeam, me, me, me, me, mmu, mmu, 0.0, 2.0 * ebeam,
			    rannums, 0, 1);
		if (dj == 0.0) {
			n_kin_fail++;
			continue;
		}

		reconstruct_momenta(&p1, &p2, &p3, &p4, &p5, &p6, &p7);
		vsub(&q1, &p1, &p3);
		vsub(&q2, &p2, &p5);
		vsub(&Q, &p6, &p7);

		m2_3 = v4dot(&p3, &p3);
		m2_5 = v4dot(&p5, &p5);
		m2_6 = v4dot(&p6, &p6);
		m2_7 = v4dot(&p7, &p7);
		if (fabs(m2_3 - me * me) > 1e-6 || fabs(m2_5 - me * me) > 1e-6 ||
		    fabs(m2_6 - mmu * mmu) > 1e-5 ||
		    fabs(m2_7 - mmu * mmu) > 1e-5) {
			n_mass_fail++;
			continue;
		}

		vadd4(&sum_in, &p1, &p2);
		vadd4(&tmp, &p3, &p5);
		vadd4(&sum_out, &tmp, &p6);
		vadd4(&sum_out, &sum_out, &p7);
		cons = fabs(sum_in.e[0] - sum_out.e[0]) +
		       fabs(sum_in.e[1] - sum_out.e[1]) +
		       fabs(sum_in.e[2] - sum_out.e[2]) +
		       fabs(sum_in.e[3] - sum_out.e[3]);
		if (cons > max_cons)
			max_cons = cons;
		if (cons > 1e-4) {
			n_cons_fail++;
			continue;
		}

		t1 = extra.t1;
		t2 = extra.t2;
		t1v = v4dot(&q1, &q1);
		t2v = v4dot(&q2, &q2);
		if (!t_match(t1, t1v) || !t_match(t2, t2v)) {
			n_t_fail++;
			continue;
		}

		per = peripp(0, 0);
		if (!isfinite(per) || per < 0.0) {
			n_peripp_bad++;
			continue;
		}

		n_ok++;

		if (soft_t(t1, t2)) {
			Bcode = civita.bb;
			sumM = eval_M11(&q1, &q2, &Q, mmu) +
			       eval_M12(&q1, &q2, &Q, &p2, mmu) +
			       eval_M21(&q1, &q2, &Q, &p1, mmu) +
			       eval_M22(&q1, &q2, &Q, &p1, &p2);
			pred = 0.25 * sumM /
			       (t1 * t2 * Bcode * t1 * t2 * Bcode);
			n_soft++;
			if (!isfinite(pred) || fabs(Bcode) < 1e-30) {
				n_soft_fail++;
			} else {
				rel = fabs(pred - per) / fmax(fabs(per), 1e-300);
				if (rel > max_soft_rel)
					max_soft_rel = rel;
				if (rel > 1e-3)
					n_soft_fail++;
			}
		}
	}

	printf("accepted=%d  kin_fail=%d  mass_fail=%d  cons_fail=%d  "
	       "t_fail=%d  peripp_bad=%d\n",
	       n_ok, n_kin_fail, n_mass_fail, n_cons_fail, n_t_fail,
	       n_peripp_bad);
	printf("max_|p_in-p_out|=%.3e\n", max_cons);
	printf("soft_|t|<1: n=%d fail=%d max_rel(peripp vs Mij)=%.3e\n", n_soft,
	       n_soft_fail, max_soft_rel);

	if (n_ok < 100) {
		fprintf(stderr, "FAIL: too few accepted points (%d)\n", n_ok);
		return 1;
	}
	if (n_mass_fail || n_cons_fail || n_t_fail || n_peripp_bad) {
		fprintf(stderr,
			"FAIL: kinematics/peripp "
			"(mass=%d cons=%d t=%d peripp=%d)\n",
			n_mass_fail, n_cons_fail, n_t_fail, n_peripp_bad);
		return 1;
	}
	if (n_soft < N_SOFT_NEED || n_soft_fail > 0) {
		fprintf(stderr,
			"FAIL: soft peripp vs ΣMij "
			"(soft=%d need=%d fail=%d)\n",
			n_soft, N_SOFT_NEED, n_soft_fail);
		return 1;
	}

	printf("PASS: check_eemm_gamgam\n");
	return 0;
}
