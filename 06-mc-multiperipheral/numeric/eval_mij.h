/*
 * Vermaseren (3.11) M_ij sum evaluator: S = (M11+M12+M21+M22) / B^2
 * with B = (q1·q2)^2 − (Q·q1)^2, Q = p6−p7.
 *
 * Conventions match scripts/mumu_mij_check.py:
 *   M12/M21 three-index εε via +eps3sq (Gram −det; opposite NPB leading minus)
 *   M22 bilinears ε(p2,q2;q1,Q) / ε(p1,q1;q2,Q)
 *
 * Hand-maintained C twin of the Python gate.
 */
#ifndef EVAL_MIJ_H
#define EVAL_MIJ_H

#include "eval_stable_35.h" /* V4, V4l, v4dot, s35_*, det helpers */

static inline double mij_eps3sq(const V4 *a, const V4 *b, const V4 *c)
{
	/* −det Gram(a,b,c); M12/M21 use +eps3sq in place of NPB (−εε) */
	return -s35_det3(v4dot(a, a), v4dot(a, b), v4dot(a, c), v4dot(b, a),
			 v4dot(b, b), v4dot(b, c), v4dot(c, a), v4dot(c, b),
			 v4dot(c, c));
}

static inline double mij_chain3(const V4 *a, const V4 *b, const V4 *c,
				const V4 *d, const V4 *e, const V4 *f)
{
	return s35_det3(v4dot(a, d), v4dot(a, e), v4dot(a, f), v4dot(b, d),
			v4dot(b, e), v4dot(b, f), v4dot(c, d), v4dot(c, e),
			v4dot(c, f));
}

static inline double eval_M11(const V4 *q1, const V4 *q2, const V4 *Q, double m)
{
	double q1q1 = v4dot(q1, q1);
	double q2q2 = v4dot(q2, q2);
	double q1q2 = v4dot(q1, q2);
	double Qq1 = v4dot(Q, q1);
	double QQ = v4dot(Q, Q);
	double B = q1q2 * q1q2 - Qq1 * Qq1;
	double mm = m * m;
	double brace =
	    B * (Qq1 * Qq1 - q1q1 * q2q2 + q1q2 * q1q2 -
		 QQ * (q1q1 + q2q2 + 2.0 * mm)) -
	    2.0 * Qq1 * Qq1 * (q1q1 + 2.0 * mm) * (q2q2 + 2.0 * mm);
	return 64.0 * q1q1 * q2q2 * brace;
}

static inline double eval_M12(const V4 *q1, const V4 *q2, const V4 *Q,
			      const V4 *p2, double m)
{
	double q1q1 = v4dot(q1, q1);
	double q1q2 = v4dot(q1, q2);
	double Qq1 = v4dot(Q, q1);
	double B = q1q2 * q1q2 - Qq1 * Qq1;
	double t1 = mij_eps3sq(p2, q2, q1) + mij_eps3sq(p2, q2, Q);
	double e_p2q2 = s35_biline(p2, q2, p2, q2);
	double e_mix = s35_biline(p2, q2, q1, Q);
	double t2 = Qq1 * Qq1 * e_p2q2 + 0.5 * e_mix * e_mix;
	return 128.0 * q1q1 * (B * t1 - (q1q1 + 2.0 * m * m) * t2);
}

static inline double eval_M21(const V4 *q1, const V4 *q2, const V4 *Q,
			      const V4 *p1, double m)
{
	double q2q2 = v4dot(q2, q2);
	double q1q2 = v4dot(q1, q2);
	double Qq1 = v4dot(Q, q1);
	double B = q1q2 * q1q2 - Qq1 * Qq1;
	double t1 = mij_eps3sq(p1, q1, q2) + mij_eps3sq(p1, q1, Q);
	double e_p1q1 = s35_biline(p1, q1, p1, q1);
	double e_mix = s35_biline(p1, q1, q2, Q);
	double t2 = Qq1 * Qq1 * e_p1q1 + 0.5 * e_mix * e_mix;
	return 128.0 * q2q2 * (B * t1 - (q2q2 + 2.0 * m * m) * t2);
}

static inline double eval_M22(const V4 *q1, const V4 *q2, const V4 *Q,
			      const V4 *p1, const V4 *p2)
{
	double q1q2 = v4dot(q1, q2);
	double Qq1 = v4dot(Q, q1);
	double B = q1q2 * q1q2 - Qq1 * Qq1;
	double Tr11 = s35_Tr(p1, q1, p1, q1);
	double Tr22 = s35_Tr(p2, q2, p2, q2);
	double Tr12 = s35_Tr(p1, q1, p2, q2);
	double EpsFull = s35_det4(p1, q1, p2, q2);
	double mid_eps = mij_chain3(p1, q1, Q, p2, q2, Q);
	double mid = mid_eps - 0.5 * Tr12 * (v4dot(Q, Q) + q1q2);
	double bil_2q1Q = s35_biline(p2, q2, q1, Q);
	double bil_1q2Q = s35_biline(p1, q1, q2, Q);
	double brace = B * (0.25 * Tr12 * Tr12 - EpsFull) - mid * mid -
		       0.125 * Tr11 * bil_2q1Q * bil_2q1Q -
		       0.125 * Tr22 * bil_1q2Q * bil_1q2Q -
		       0.25 * Qq1 * Qq1 * Tr11 * Tr22;
	return 512.0 * brace;
}

static inline double eval_mij_sum_over_B2(const V4 *q1, const V4 *q2,
					 const V4 *Q, const V4 *p1,
					 const V4 *p2, double m)
{
	double q1q2 = v4dot(q1, q2);
	double Qq1 = v4dot(Q, q1);
	double B = q1q2 * q1q2 - Qq1 * Qq1;
	double sum = eval_M11(q1, q2, Q, m) + eval_M12(q1, q2, Q, p2, m) +
		     eval_M21(q1, q2, Q, p1, m) + eval_M22(q1, q2, Q, p1, p2);
	return sum / (B * B);
}

/* ---- long double ---- */

static inline long double mij_eps3sq_l(const V4l *a, const V4l *b, const V4l *c)
{
	return -s35_det3_l(v4dot_l(a, a), v4dot_l(a, b), v4dot_l(a, c),
			   v4dot_l(b, a), v4dot_l(b, b), v4dot_l(b, c),
			   v4dot_l(c, a), v4dot_l(c, b), v4dot_l(c, c));
}

static inline long double mij_chain3_l(const V4l *a, const V4l *b, const V4l *c,
					const V4l *d, const V4l *e, const V4l *f)
{
	return s35_det3_l(v4dot_l(a, d), v4dot_l(a, e), v4dot_l(a, f),
			  v4dot_l(b, d), v4dot_l(b, e), v4dot_l(b, f),
			  v4dot_l(c, d), v4dot_l(c, e), v4dot_l(c, f));
}

static inline long double eval_M11_l(const V4l *q1, const V4l *q2, const V4l *Q,
				     long double m)
{
	long double q1q1 = v4dot_l(q1, q1);
	long double q2q2 = v4dot_l(q2, q2);
	long double q1q2 = v4dot_l(q1, q2);
	long double Qq1 = v4dot_l(Q, q1);
	long double QQ = v4dot_l(Q, Q);
	long double B = q1q2 * q1q2 - Qq1 * Qq1;
	long double mm = m * m;
	long double brace =
	    B * (Qq1 * Qq1 - q1q1 * q2q2 + q1q2 * q1q2 -
		 QQ * (q1q1 + q2q2 + 2.0L * mm)) -
	    2.0L * Qq1 * Qq1 * (q1q1 + 2.0L * mm) * (q2q2 + 2.0L * mm);
	return 64.0L * q1q1 * q2q2 * brace;
}

static inline long double eval_M12_l(const V4l *q1, const V4l *q2, const V4l *Q,
				     const V4l *p2, long double m)
{
	long double q1q1 = v4dot_l(q1, q1);
	long double q1q2 = v4dot_l(q1, q2);
	long double Qq1 = v4dot_l(Q, q1);
	long double B = q1q2 * q1q2 - Qq1 * Qq1;
	long double t1 = mij_eps3sq_l(p2, q2, q1) + mij_eps3sq_l(p2, q2, Q);
	long double e_p2q2 = s35_biline_l(p2, q2, p2, q2);
	long double e_mix = s35_biline_l(p2, q2, q1, Q);
	long double t2 = Qq1 * Qq1 * e_p2q2 + 0.5L * e_mix * e_mix;
	return 128.0L * q1q1 * (B * t1 - (q1q1 + 2.0L * m * m) * t2);
}

static inline long double eval_M21_l(const V4l *q1, const V4l *q2, const V4l *Q,
				     const V4l *p1, long double m)
{
	long double q2q2 = v4dot_l(q2, q2);
	long double q1q2 = v4dot_l(q1, q2);
	long double Qq1 = v4dot_l(Q, q1);
	long double B = q1q2 * q1q2 - Qq1 * Qq1;
	long double t1 = mij_eps3sq_l(p1, q1, q2) + mij_eps3sq_l(p1, q1, Q);
	long double e_p1q1 = s35_biline_l(p1, q1, p1, q1);
	long double e_mix = s35_biline_l(p1, q1, q2, Q);
	long double t2 = Qq1 * Qq1 * e_p1q1 + 0.5L * e_mix * e_mix;
	return 128.0L * q2q2 * (B * t1 - (q2q2 + 2.0L * m * m) * t2);
}

static inline long double eval_M22_l(const V4l *q1, const V4l *q2, const V4l *Q,
				     const V4l *p1, const V4l *p2)
{
	long double q1q2 = v4dot_l(q1, q2);
	long double Qq1 = v4dot_l(Q, q1);
	long double B = q1q2 * q1q2 - Qq1 * Qq1;
	long double Tr11 = s35_Tr_l(p1, q1, p1, q1);
	long double Tr22 = s35_Tr_l(p2, q2, p2, q2);
	long double Tr12 = s35_Tr_l(p1, q1, p2, q2);
	long double EpsFull = s35_det4_l(p1, q1, p2, q2);
	long double mid_eps = mij_chain3_l(p1, q1, Q, p2, q2, Q);
	long double mid = mid_eps - 0.5L * Tr12 * (v4dot_l(Q, Q) + q1q2);
	long double bil_2q1Q = s35_biline_l(p2, q2, q1, Q);
	long double bil_1q2Q = s35_biline_l(p1, q1, q2, Q);
	long double brace = B * (0.25L * Tr12 * Tr12 - EpsFull) - mid * mid -
			    0.125L * Tr11 * bil_2q1Q * bil_2q1Q -
			    0.125L * Tr22 * bil_1q2Q * bil_1q2Q -
			    0.25L * Qq1 * Qq1 * Tr11 * Tr22;
	return 512.0L * brace;
}

static inline long double eval_mij_sum_over_B2_l(const V4l *q1, const V4l *q2,
						 const V4l *Q, const V4l *p1,
						 const V4l *p2, long double m)
{
	long double q1q2 = v4dot_l(q1, q2);
	long double Qq1 = v4dot_l(Q, q1);
	long double B = q1q2 * q1q2 - Qq1 * Qq1;
	long double sum = eval_M11_l(q1, q2, Q, m) + eval_M12_l(q1, q2, Q, p2, m) +
			  eval_M21_l(q1, q2, Q, p1, m) +
			  eval_M22_l(q1, q2, Q, p1, p2);
	return sum / (B * B);
}

#endif /* EVAL_MIJ_H */
