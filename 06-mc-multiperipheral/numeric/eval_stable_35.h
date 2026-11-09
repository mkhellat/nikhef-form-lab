/*
 * Stable Vermaseren (3.5) evaluator for γ*γ*→μ⁺μ⁻ (subreaction).
 * Matches scripts/mumu_stable_35_check.py (incl. minus on last bilinears).
 * Q = p6 − p7. Minkowski (+−−−).
 *
 * Hand-maintained C twin of the Python gate; keep in sync if (3.5) signs
 * change. Double and long-double entry points for digit-loss scans.
 */
#ifndef EVAL_STABLE_35_H
#define EVAL_STABLE_35_H

#include <math.h>

typedef struct {
	double e[4];
} V4;

static inline double v4dot(const V4 *a, const V4 *b)
{
	return a->e[0] * b->e[0] - a->e[1] * b->e[1] - a->e[2] * b->e[2] -
	       a->e[3] * b->e[3];
}

static inline double s35_Tr(const V4 *qi, const V4 *ei, const V4 *qj,
			    const V4 *ej)
{
	return 2.0 * (v4dot(qi, qj) * v4dot(ei, ej) - v4dot(qi, ej) * v4dot(ei, qj));
}

static inline double s35_biline(const V4 *a, const V4 *b, const V4 *c,
				const V4 *d)
{
	return -2.0 * (v4dot(a, c) * v4dot(b, d) - v4dot(a, d) * v4dot(b, c));
}

static inline double s35_det3(double m00, double m01, double m02, double m10,
			      double m11, double m12, double m20, double m21,
			      double m22)
{
	return m00 * (m11 * m22 - m12 * m21) - m01 * (m10 * m22 - m12 * m20) +
	       m02 * (m10 * m21 - m11 * m20);
}

static inline double s35_det4(const V4 *a, const V4 *b, const V4 *c,
			      const V4 *d)
{
	/* det Gram(a,b,c,d) = ε^{abcd}ε_{abcd} in (+−−−) as used in Python */
	double G[4][4];
	const V4 *v[4] = {a, b, c, d};
	int i, j;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			G[i][j] = v4dot(v[i], v[j]);
	/* Laplace expansion on first row */
	return G[0][0] * s35_det3(G[1][1], G[1][2], G[1][3], G[2][1], G[2][2],
				  G[2][3], G[3][1], G[3][2], G[3][3]) -
	       G[0][1] * s35_det3(G[1][0], G[1][2], G[1][3], G[2][0], G[2][2],
				  G[2][3], G[3][0], G[3][2], G[3][3]) +
	       G[0][2] * s35_det3(G[1][0], G[1][1], G[1][3], G[2][0], G[2][1],
				  G[2][3], G[3][0], G[3][1], G[3][3]) -
	       G[0][3] * s35_det3(G[1][0], G[1][1], G[1][2], G[2][0], G[2][1],
				  G[2][2], G[3][0], G[3][1], G[3][2]);
}

static inline double s35_Q12Q(const V4 *q1, const V4 *e1, const V4 *Q,
			      const V4 *q2, const V4 *e2)
{
	/* ε^{q1 e1 Q μ} ε_{q2 e2 μ Q} via 3×3 Gram (Python Q12Q) */
	return s35_det3(v4dot(q1, q2), v4dot(q1, e2), v4dot(q1, Q),
			v4dot(e1, q2), v4dot(e1, e2), v4dot(e1, Q),
			v4dot(Q, q2), v4dot(Q, e2), v4dot(Q, Q));
}

/* S_(3.5) = Num / B^2 with B = (q1·q2)^2 − (Q·q1)^2 */
static inline double eval_stable_35(const V4 *q1, const V4 *q2, const V4 *e1,
				    const V4 *e2, const V4 *Q)
{
	double q1q2 = v4dot(q1, q2);
	double Qq1 = v4dot(Q, q1);
	double B = q1q2 * q1q2 - Qq1 * Qq1;
	double Tr11 = s35_Tr(q1, e1, q1, e1);
	double Tr22 = s35_Tr(q2, e2, q2, e2);
	double Tr12 = s35_Tr(q1, e1, q2, e2);
	double EpsFull = s35_det4(q1, e1, q2, e2);
	double q12q = s35_Q12Q(q1, e1, Q, q2, e2);
	double Eps2Q = s35_biline(q2, e2, q1, Q);
	double Eps1Q = s35_biline(q1, e1, q2, Q);
	double mid = q12q - 0.5 * Tr12 * (v4dot(Q, Q) + q1q2);
	double Num = 8.0 * B * (0.25 * Tr12 * Tr12 - EpsFull) -
		     2.0 * Qq1 * Qq1 * Tr11 * Tr22 - 8.0 * mid * mid -
		     Tr11 * Eps2Q * Eps2Q - Tr22 * Eps1Q * Eps1Q;
	return Num / (B * B);
}

/* ---- long double twin ---- */

typedef struct {
	long double e[4];
} V4l;

static inline long double v4dot_l(const V4l *a, const V4l *b)
{
	return a->e[0] * b->e[0] - a->e[1] * b->e[1] - a->e[2] * b->e[2] -
	       a->e[3] * b->e[3];
}

static inline long double s35_Tr_l(const V4l *qi, const V4l *ei, const V4l *qj,
				   const V4l *ej)
{
	return 2.0L * (v4dot_l(qi, qj) * v4dot_l(ei, ej) -
		       v4dot_l(qi, ej) * v4dot_l(ei, qj));
}

static inline long double s35_biline_l(const V4l *a, const V4l *b, const V4l *c,
				       const V4l *d)
{
	return -2.0L * (v4dot_l(a, c) * v4dot_l(b, d) - v4dot_l(a, d) * v4dot_l(b, c));
}

static inline long double s35_det3_l(long double m00, long double m01,
				     long double m02, long double m10,
				     long double m11, long double m12,
				     long double m20, long double m21,
				     long double m22)
{
	return m00 * (m11 * m22 - m12 * m21) - m01 * (m10 * m22 - m12 * m20) +
	       m02 * (m10 * m21 - m11 * m20);
}

static inline long double s35_det4_l(const V4l *a, const V4l *b, const V4l *c,
				     const V4l *d)
{
	long double G[4][4];
	const V4l *v[4] = {a, b, c, d};
	int i, j;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			G[i][j] = v4dot_l(v[i], v[j]);
	return G[0][0] * s35_det3_l(G[1][1], G[1][2], G[1][3], G[2][1], G[2][2],
				    G[2][3], G[3][1], G[3][2], G[3][3]) -
	       G[0][1] * s35_det3_l(G[1][0], G[1][2], G[1][3], G[2][0], G[2][2],
				    G[2][3], G[3][0], G[3][2], G[3][3]) +
	       G[0][2] * s35_det3_l(G[1][0], G[1][1], G[1][3], G[2][0], G[2][1],
				    G[2][3], G[3][0], G[3][1], G[3][3]) -
	       G[0][3] * s35_det3_l(G[1][0], G[1][1], G[1][2], G[2][0], G[2][1],
				    G[2][2], G[3][0], G[3][1], G[3][2]);
}

static inline long double s35_Q12Q_l(const V4l *q1, const V4l *e1, const V4l *Q,
				     const V4l *q2, const V4l *e2)
{
	return s35_det3_l(v4dot_l(q1, q2), v4dot_l(q1, e2), v4dot_l(q1, Q),
			  v4dot_l(e1, q2), v4dot_l(e1, e2), v4dot_l(e1, Q),
			  v4dot_l(Q, q2), v4dot_l(Q, e2), v4dot_l(Q, Q));
}

static inline long double eval_stable_35_l(const V4l *q1, const V4l *q2,
					   const V4l *e1, const V4l *e2,
					   const V4l *Q)
{
	long double q1q2 = v4dot_l(q1, q2);
	long double Qq1 = v4dot_l(Q, q1);
	long double B = q1q2 * q1q2 - Qq1 * Qq1;
	long double Tr11 = s35_Tr_l(q1, e1, q1, e1);
	long double Tr22 = s35_Tr_l(q2, e2, q2, e2);
	long double Tr12 = s35_Tr_l(q1, e1, q2, e2);
	long double EpsFull = s35_det4_l(q1, e1, q2, e2);
	long double q12q = s35_Q12Q_l(q1, e1, Q, q2, e2);
	long double Eps2Q = s35_biline_l(q2, e2, q1, Q);
	long double Eps1Q = s35_biline_l(q1, e1, q2, Q);
	long double mid = q12q - 0.5L * Tr12 * (v4dot_l(Q, Q) + q1q2);
	long double Num = 8.0L * B * (0.25L * Tr12 * Tr12 - EpsFull) -
			  2.0L * Qq1 * Qq1 * Tr11 * Tr22 - 8.0L * mid * mid -
			  Tr11 * Eps2Q * Eps2Q - Tr22 * Eps1Q * Eps1Q;
	return Num / (B * B);
}

static inline void v4_to_l(V4l *o, const V4 *a)
{
	o->e[0] = (long double)a->e[0];
	o->e[1] = (long double)a->e[1];
	o->e[2] = (long double)a->e[2];
	o->e[3] = (long double)a->e[3];
}

#endif /* EVAL_STABLE_35_H */
