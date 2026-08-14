/*
	Numerical confirmation of the s-scaling claim in uam19 part5.pdf's
	discussion of formula (3.2) (Vermaseren 1983, p.353): "the first term
	combines with 1/(t1 t2)^2 which by itself would cause problems with
	physics because it would give a crosssection that grows with powers
	of s."

	This driver instruments the REAL kinc2 pickin.c/orient.c/pi0.c code
	directly -- it is not a reimplementation of the kinematics or the
	matrix element. It scans the beam energy sqrt(s) while holding the
	photon virtualities t1, t2 FIXED at a small (near-real-photon) target
	value via bisection on the random numbers that pickin.c's mapt()
	uses to select them, then reads out pi0.c's four Levi-Civita/Gram
	terms (part1..part4, matching formula (3.2)'s four terms exactly)
	and their sum.

	Result (see ../notes/topic5-s-scaling-pi0.md for the full discussion):
	part1 -- the term derived analytically in scripts/gram_s_scaling.py
	to scale like s^2/(t1 t2) at fixed t1, t2 -- dominates at every s
	tested, and the SUM tracks part1 with no cancellation. This confirms,
	against the actual code, that the four terms of formula (3.2) do NOT
	cancel each other down to a decent (non-growing) s-behavior; the
	paper says so explicitly (p.353: "all 4 terms are positive so no
	cancellations occur"). The "growth with powers of s would be bad"
	warning in part5.pdf is about a DIFFERENT, hypothetical computation:
	if you contracted the Levi-Civita tensors into ordinary 4-vector dot
	products (rather than using pickin.c's stable factored form,
	`levi.gram = (1-yy4^2)*dd/ap`), you would get individual dot-product
	terms of size O(s^2 Mpi^4) that must cancel against EACH OTHER
	(within that single dot-product expansion, not across part1..part4)
	down to something O(t1 t2) in size -- a numerically catastrophic
	cancellation that pickin.c avoids entirely by never performing that
	expansion.

	Depends on the vendored kinc2 sources in ../kinc2/ (pickin.c, orient.c,
	mapt.c, mapw.c, mapla.c, ipow.c, iipow.c, ranf.c, pi0.c, axolib.h --
	copied verbatim from the kinc2.tar.gz archive referenced throughout
	this directory's notes; see ../kinc2/README.md. This program is
	buildable and runnable directly from the repo, with no external
	paths required, but is not itself part of that vendored code.)

	Build (from this directory, kinc2-driver/):
		cc -Wall -g -I../kinc2 -c ../kinc2/orient.c ../kinc2/pickin.c \
			../kinc2/mapt.c ../kinc2/mapw.c ../kinc2/mapla.c \
			../kinc2/ipow.c ../kinc2/iipow.c ../kinc2/ranf.c \
			../kinc2/pi0.c
		cc -Wall -g -I../kinc2 -c scan_pi0_scaling.c
		cc -o scan_pi0_scaling *.o -lm
		./scan_pi0_scaling
*/

#include "axolib.h"

double pi0(int par);

/* pi0.c reads eemminput.mu (the pi0 mass) as an extern global, normally
   defined and set up in eee.c's main(). We provide the storage here
   directly instead of linking eee.o, which pulls in the VEGAS/plotting
   machinery we don't need for this scan. */
EEMMINPUT eemminput;

static double me = 0.000511;   /* electron mass, GeV */
static double mpi = 0.135;     /* pi0 mass, GeV */

/* One kinematics + matrix-element evaluation at the given s and random
   numbers. Returns the pickin/orient Jacobian (0 signals a kinematically
   forbidden point for this s and rannums). */
static double try_point(double s, double *rannums)
{
	return orient(s, me, me, me, mpi, me, rannums, 0);
}

/*
	Bisects rannums[idx] in (0,1) so that the resulting extra.t1 (if
	which==1) or extra.t2 (if which==2) equals `target`.

	Why bisection and not a closed-form inversion: pickin.c's mapt()
	maps rannums[0] -> t1 (and rannums[1] -> t2) via a monotonic
	log-interpolation t(x) = tmin*(tmax/tmin)^x between s-dependent
	bounds t1min(s), t1max(s) computed deep inside pickin() (pickin.c
	lines 90-92, coupled to t2's own bounds via s2, s1, and several
	Gram-related intermediate quantities). Reproducing that formula
	independently would risk exactly the kind of "derive kinematics by
	hand instead of trusting the verified code" mistake this whole
	exercise is checking for. Bisecting against the real orient() output
	is monotonic (mapt is monotonic in x) and non-invasive: it finds the
	x that gives the target t without assuming anything about the bounds'
	functional form.
*/
static int bisect_for_target(double s, double *rannums, int idx, int which, double target)
{
	double lo = 1e-6, hi = 1 - 1e-6, mid = 0.5;
	double dj, t_lo, t_hi, tval;
	int it;

	rannums[idx] = lo;
	dj = try_point(s, rannums);
	if (dj == 0) return 0;
	t_lo = (which == 1) ? extra.t1 : extra.t2;

	rannums[idx] = hi;
	dj = try_point(s, rannums);
	if (dj == 0) return 0;
	t_hi = (which == 1) ? extra.t1 : extra.t2;

	if (target < (t_lo < t_hi ? t_lo : t_hi) || target > (t_lo > t_hi ? t_lo : t_hi)) {
		return 0;  /* target t not achievable in this s's [t_lo, t_hi] range */
	}

	for (it = 0; it < 80; it++) {
		mid = 0.5 * (lo + hi);
		rannums[idx] = mid;
		dj = try_point(s, rannums);
		if (dj == 0) {
			/* Kinematics failed exactly at this trial point (can happen
			   near a boundary); nudge the search inward and retry. */
			lo = (t_lo < t_hi) ? mid : lo;
			continue;
		}
		tval = (which == 1) ? extra.t1 : extra.t2;
		if ((t_lo < t_hi && tval < target) || (t_lo > t_hi && tval > target)) {
			lo = mid; t_lo = tval;
		} else {
			hi = mid; t_hi = tval;
		}
		if (fabs(tval - target) < 1e-9 * fabs(target)) break;
	}
	rannums[idx] = mid;
	return 1;
}

int main(void)
{
	double rannums[4];
	/* sqrt(s) values spanning three orders of magnitude in s, all with
	   the same fixed small photon virtualities below. */
	double sqrts_vals[] = {0.5, 1.0, 3.0, 10.0, 30.0, 3000.0, 10000.0};
	double t1_target = -0.001;  /* fixed small spacelike virtuality, GeV^2 */
	double t2_target = -0.001;
	int i;

	eemminput.mu = mpi;

	printf("%12s %14s %14s %14s %14s %14s %14s %14s %14s\n",
		"sqrt(s)", "s", "t1", "t2", "part1", "part2", "part3", "part4", "sum");

	for (i = 0; i < (int)(sizeof(sqrts_vals) / sizeof(sqrts_vals[0])); i++) {
		double sqrts = sqrts_vals[i];
		double s = sqrts * sqrts;
		double dj, part1, part2, part3, part4, sum, tt, la;

		rannums[2] = 0.5;  /* s2 selection: mid-range, not the focus of this scan */
		rannums[3] = 0.5;  /* azimuthal cos(yy4): mid-range */
		rannums[0] = 0.5;  /* initial guess before bisection */
		rannums[1] = 0.5;

		if (!bisect_for_target(s, rannums, 0, 1, t1_target)) {
			printf("%12.4g  (t1 target unreachable at this s)\n", sqrts);
			continue;
		}
		if (!bisect_for_target(s, rannums, 1, 2, t2_target)) {
			printf("%12.4g  (t2 target unreachable at this s)\n", sqrts);
			continue;
		}

		dj = try_point(s, rannums);  /* final evaluation at the converged rannums */
		if (dj == 0) {
			printf("%12.4g  (kinematics failed at converged point)\n", sqrts);
			continue;
		}

		/* This block is pi0.c's pi0() function, inlined so we can print
		   the individual parts (pi0() itself only returns the sum). */
		tt = extra.t1 * extra.t2;
		la = (mpi * mpi - extra.t1 - extra.t2) * 0.5;
		la = -(la * la - tt);
		part1 = -64 * levi.gram / (tt * tt);
		part2 = -16 * extra.t1 * levi.dd2 / (tt * tt);
		part3 = -16 * extra.t2 * levi.dd4 / (tt * tt);
		part4 = -4 * tt * la / (tt * tt);
		sum = part1 + part2 + part3 + part4;

		printf("%12.4g %14.6g %14.6g %14.6g %14.6g %14.6g %14.6g %14.6g %14.6g\n",
			sqrts, s, extra.t1, extra.t2, part1, part2, part3, part4, sum);
	}
	return 0;
}
