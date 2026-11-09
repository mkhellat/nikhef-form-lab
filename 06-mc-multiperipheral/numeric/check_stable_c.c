/*
 * Smoke gate: C eval_stable_35 and eval_mij match fixed Python reference
 * values (same kinematics as scripts/mumu_*_check.py conventions).
 *
 * Build: make check && ./check_stable_c
 */

#include "eval_mij.h"
#include "eval_stable_35.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int check_rel(const char *name, double got, double want, double tol)
{
	double rel = fabs(got - want) / fmax(fabs(want), 1e-300);
	if (rel >= tol) {
		fprintf(stderr, "FAIL %s: got=%.12e want=%.12e rel=%.3e\n", name,
			got, want, rel);
		return 1;
	}
	printf("ok %s rel=%.3e\n", name, rel);
	return 0;
}

int main(void)
{
	/* Fixed kinematics: W=20, eps=1e-2, delta=0.1, m=muon, e2=yhat */
	const double W = 20.0;
	const double m = 0.105658;
	const double eps = 1e-2;
	const double delta = 1e-1;
	double E = W / 2.0;
	double qz = sqrt(E * E + eps);
	double p = sqrt(E * E - m * m);
	V4 q1 = {{E, 0, 0, qz}};
	V4 q2 = {{E, 0, 0, -qz}};
	V4 p6 = {{E, p, 0, 0}};
	V4 p7 = {{E, -p, 0, 0}};
	V4 Q = {{p6.e[0] - p7.e[0], p6.e[1] - p7.e[1], p6.e[2] - p7.e[2],
		 p6.e[3] - p7.e[3]}};
	V4 e1 = {{q1.e[0], q1.e[1] + delta, q1.e[2], q1.e[3]}};
	V4 e2 = {{0, 0, 1, 0}};
	/* massless electrons with 2 p·q = q² along soft transverse */
	V4 p1 = {{q1.e[0] * 0.5, 0.3, 0.0, q1.e[3] * 0.5}};
	V4 p2 = {{q2.e[0] * 0.5, 0.0, 0.3, q2.e[3] * 0.5}};
	int rc = 0;

	/* Python S_35 at this point (verified 2026-11-09) */
	rc |= check_rel("S_35", eval_stable_35(&q1, &q2, &e1, &e2, &Q),
			0.07999200149300065, 1e-12);

	/* Same p1,p2 as above: Python ΣM_ij/B^2 (verified 2026-11-09) */
	{
		V4l q1l, q2l, Ql, p1l, p2l;
		double Sd;
		long double Sl;
		v4_to_l(&q1l, &q1);
		v4_to_l(&q2l, &q2);
		v4_to_l(&Ql, &Q);
		v4_to_l(&p1l, &p1);
		v4_to_l(&p2l, &p2);
		Sd = eval_mij_sum_over_B2(&q1, &q2, &Q, &p1, &p2, m);
		Sl = eval_mij_sum_over_B2_l(&q1l, &q2l, &Ql, &p1l, &p2l,
					    (long double)m);
		rc |= check_rel("Mij_vs_py", Sd, 4.383602384069826, 1e-12);
		rc |= check_rel("Mij_d_vs_ld", Sd, (double)Sl, 1e-10);
	}

	if (rc)
		return 1;
	printf("PASS: check_stable_c\n");
	return 0;
}
