/*
	VEGAS scan of the integrated e-e+ -> e-e+ pi0 cross-section proxy
	versus sqrt(s), with optional |t| and scattering-angle cuts, compared
	with the leading-log (LL) phase-space volume law.

	WHY THIS DRIVER EXISTS.  mc_xsec_s_scan.c integrates the same
	integrand with plain random points on the unit cube.  That is
	already importance sampled: pickin() maps both t's with mapt(...,1),
	i.e. dt/t (../kinc2/pickin.c:93,142), which flattens most of the
	1/(t1 t2) peak, and its 2M-point numbers agree with this driver to
	about 0.1 %.  What VEGAS adds is an error estimate and a chi^2 per
	run, a further adaptive gain (~1e-4 relative error with 2e5 points
	per iteration), and efficiency under cuts that veto most of the
	dt/t-sampled points (--thmin).  That lets it
	resolve the *shape* of sigma(s): that is what the notes (see
	../notes/topic5-s-scaling-pi0.md) need in order to separate

	  - the uncut, inclusive growth   (a polynomial in ln s),
	  - the growth with a fixed |t| >= t0 cut  (a weaker polynomial), and
	  - the s-independent plateau once the electrons are kept away from
	    the beam by an angular cut (--thmin).

	INTEGRAND (identical conventions to mc_xsec_s_scan.c, no
	re-derivation of any kinematics or matrix element; orient(), pickin()
	and pi0() are the vendored ../kinc2/ routines):

	    f(r) = dj(r) * |M|^2_pi0 / (2 s)         r in the 4D unit cube

	where dj is the phase-space Jacobian returned by orient() and
	|M|^2 comes from pi0(0).  The flux factor 1/(2s) is the one in
	eee.c's overallconstant.  The coupling uses a FIXED alpha (nothing
	in this integrand runs alpha), unlike eee.c, so the s-dependence
	below is pure phase-space/matrix-element and not vacuum polarisation.
	The output is called sigma_raw and is in "raw/2s" units.

	NO ABSOLUTE CROSS SECTION IS PRINTED.  Two normalisation issues are
	still open and are deliberately not papered over here:
	  (a) an unresolved overall factor 2 between the coupling and the
	      spin-sum conventions for the pi0 matrix element, and
	  (b) a missing (2 pi)^-5 phase-space factor in eee.c's
	      overallconstant.
	Until both are settled, only the *shape* in s and the ratio to the
	LL law below are meaningful.  Do not convert sigma_raw to pb.

	CUTS (applied as a veto right after orient(), so cut points simply
	contribute zero; vegas then adapts to the allowed region):
	  --tmin T    both |t1| and |t2| >= T  (GeV^2)
	  --tmax T    both |t1| and |t2| <= T  (GeV^2)
	  --thmin TH  both scattered leptons have polar angle >= TH (rad)
	              w.r.t. their own beam direction.  The electron angle is
	              atan2(st3, ct3); the positron one is atan2(st5, -ct5)
	              because the positron travels along -z.
	The default is no cut at all.

	LEADING-LOG LAW.  Write M for the pi0 mass, m_e for the electron
	mass, A = ln(M^2/m_e^2) and L = ln(s/M^2).  Each photon propagator
	1/t^2 integrated over its own virtuality gives a logarithm
	(this is the Weizsaecker-Williams/equivalent-photon picture, cf.
	Budnev, Ginzburg, Meledin, Serbo, Phys. Rep. 15 (1975) 181; and
	Vermaseren, Nucl. Phys. B229 (1983) 347 for the 2->3 form used by
	kinc2).  Counting the volume of the (ln|t1|, ln|t2|, ln s_gg)
	region that the integrand covers uniformly gives

	  uncut                V = A^2 L + 3 A L^2 + (5/3) L^3
	  |t| <= M^2           V = A^2 L + 2 A L^2 + (2/3) L^3
	  |t| >= t0 (fixed)    V = B^2 L +   B L^2 +   L^3/3 , B = ln(M^2/t0)

	The printed column  ratio = sigma_raw / V_LL  is NOT divided by any
	constant C: the prototype fits do not define one (it would absorb the
	unresolved normalisation above).  So the *level* of the ratio means
	nothing; only its s-dependence does.  The checks below use only that
	dependence (ratio must increase with s in the uncut scan).
	For --thmin (and any cut combination not listed above) there is no
	closed form we trust, so V_LL and ratio are printed as nan.  The
	angular cut is instead characterised by sigma_raw being flat in s
	(the printed sigma column itself): the cut removes the collinear
	logarithms so the leading ln(s) growth is gone.

	CHECK MODE (--check; used by `make check-xsec-pi0`, exit status 0/1):
	  G5  --thmin 0.030 at sqrt(s) = 30, 100, 1000 GeV: sigma flat, each
	      within max(3 combined sd, 0.5 %) of the 1000 GeV value
	      (long runs give 1300.8 +- 0.4 from 30 GeV to 100 TeV).
	  G2  uncut sigma at 1, 10, 100 GeV within 0.3 % of the
	      plain-random-point (dt/t-mapped) numbers from
	      mc_xsec_s_scan (2M points, srand 12345, run once and
	      pinned below).
	  G4  uncut sigma increasing with s, and sigma/V_LL increasing,
	      over the default grid 1..1000 GeV, 7 points.

	Options: --lo --hi --n  (log-spaced sqrt(s) grid in GeV, default
	1, 1000, 7), --calls N (per iteration, default 200000), --iters N
	(default 10), --seed N (seeds the vendored ranf via rcargo; default
	is ranf's own fixed seed), --check.

	References for method: G. P. Lepage, J. Comput. Phys. 27 (1978) 192
	(VEGAS, vendored in ../vegas/vegas.c).

	Build (from this directory, kinc2-driver/; or just `make
	xsec_pi0_vegas`):
		cc -O2 -Wall -I../kinc2 -I../vegas -o xsec_pi0_vegas \
			xsec_pi0_vegas.c ../vegas/vegas.c \
			../kinc2/orient.c ../kinc2/pickin.c ../kinc2/mapt.c \
			../kinc2/mapw.c ../kinc2/mapla.c ../kinc2/ipow.c \
			../kinc2/iipow.c ../kinc2/ranf.c ../kinc2/pi0.c -lm
		./xsec_pi0_vegas --thmin 0.030
*/

#include "axolib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double pi0(int par);

EEMMINPUT eemminput;

static const double me = 0.000511;	/* electron mass, GeV */
static const double mpi = 0.135;	/* pi0 mass, GeV */

/* Cuts live in file-scope statics because vegas() calls fcn(double *)
   with no user-data pointer.  Zero means "cut off". */
static double s_now;			/* s for the current point */
static double cut_tmin, cut_tmax, cut_thmin;

static double fcn(double *r)
{
	double dj = orient(s_now, me, me, me, mpi, me, r, 0);
	double m2;

	if (dj == 0)
		return 0;		/* kinematically forbidden */
	if (cut_tmin > 0 && (-extra.t1 < cut_tmin || -extra.t2 < cut_tmin))
		return 0;
	if (cut_tmax > 0 && (-extra.t1 > cut_tmax || -extra.t2 > cut_tmax))
		return 0;
	if (cut_thmin > 0) {
		/* orient() leaves the sines/cosines in `variab`.  The
		   positron (5) goes along -z, hence the minus sign. */
		double th3 = atan2(variab.st3, variab.ct3);
		double th5 = atan2(variab.st5, -variab.ct5);

		if (th3 < cut_thmin || th5 < cut_thmin)
			return 0;
	}
	m2 = pi0(0);
	if (m2 < 0)
		return 0;		/* eepi.c also drops these */
	return dj * m2 / (2 * s_now);
}

/* One VEGAS run at a given sqrt(s). */
static void run_point(double sq, long calls, int iters,
		      double *sig, double *err, double *chi2)
{
	s_now = sq * sq;
	/* accuracy 1e-5 is never reached, so all iters always run */
	vegas(fcn, 1e-5, 4, calls, iters, 0);
	*sig = avgi;
	*err = sd;
	*chi2 = chi2a;
}

/* Leading-log volume; returns nan when no closed form applies. */
static double vll(double s)
{
	double A = log(mpi * mpi / (me * me));
	double L = log(s / (mpi * mpi));

	if (cut_thmin > 0)
		return NAN;
	if (cut_tmin == 0 && cut_tmax == 0)
		return A * A * L + 3 * A * L * L + 5.0 / 3.0 * L * L * L;
	if (cut_tmin == 0 && fabs(cut_tmax / (mpi * mpi) - 1) < 1e-9)
		return A * A * L + 2 * A * L * L + 2.0 / 3.0 * L * L * L;
	if (cut_tmax == 0) {
		double B = log(mpi * mpi / cut_tmin);

		return B * B * L + B * L * L + L * L * L / 3;
	}
	return NAN;
}

static double grid_point(double lo, double hi, int n, int i)
{
	if (n < 2)
		return lo;
	return lo * pow(hi / lo, (double)i / (n - 1));
}

/* --check helpers */
static int gate(const char *name, int ok)
{
	printf("  %s %s\n", ok ? "ok  " : "BAD ", name);
	return ok;
}

static int do_check(long calls, int iters)
{
	/* Reference: mc_xsec_s_scan (plain random points, dt/t-mapped
	   by pickin, 2M points, srand(12345)),
	   column raw_avg/(2s), run on this worktree. */
	static const double ref_sq[3] = {1.0, 10.0, 100.0};
	static const double ref[3] = {5.467718e+03, 3.172604e+04,
				      8.761212e+04};
	static const double g5_sq[3] = {30.0, 100.0, 1000.0};
	double sig[7], err[7], chi2, sg[3], eg[3];
	double lo = 1, hi = 1000, ratio[7];
	int i, ok = 1, good;

	/* G5: angular cut -> flat */
	cut_tmin = cut_tmax = 0;
	cut_thmin = 0.030;
	for (i = 0; i < 3; i++)
		run_point(g5_sq[i], calls, iters, &sg[i], &eg[i], &chi2);
	for (i = 0; i < 2; i++) {
		double d = fabs(sg[i] - sg[2]);
		double tol = 3 * sqrt(eg[i] * eg[i] + eg[2] * eg[2]);
		char buf[128];

		if (tol < 0.005 * sg[2])
			tol = 0.005 * sg[2];
		snprintf(buf, sizeof buf,
			 "G5 thmin=0.030 sqrt_s=%g: %.2f vs %.2f (|d|=%.2f tol=%.2f)",
			 g5_sq[i], sg[i], sg[2], d, tol);
		ok &= gate(buf, d <= tol);
	}

	/* G2 and G4: uncut default grid */
	cut_thmin = 0;
	for (i = 0; i < 7; i++) {
		double sq = grid_point(lo, hi, 7, i);

		run_point(sq, calls, iters, &sig[i], &err[i], &chi2);
		ratio[i] = sig[i] / vll(sq * sq);
	}
	for (i = 0; i < 3; i++) {
		int k = (i == 0) ? 0 : (i == 1 ? 2 : 4);  /* 1,10,100 GeV */
		double d = fabs(sig[k] / ref[i] - 1);
		char buf[128];

		snprintf(buf, sizeof buf,
			 "G2 uncut sqrt_s=%g: %.1f vs ref %.1f (%.3f%%)",
			 ref_sq[i], sig[k], ref[i], 100 * d);
		ok &= gate(buf, d <= 0.003);
	}
	good = 1;
	for (i = 1; i < 7; i++)
		if (!(sig[i] > sig[i - 1]))
			good = 0;
	ok &= gate("G4 uncut sigma increases with s", good);
	good = 1;
	for (i = 1; i < 7; i++)
		if (!(ratio[i] > ratio[i - 1]))
			good = 0;
	ok &= gate("G4 uncut sigma/V_LL increases with s", good);

	printf("%s: xsec_pi0_vegas\n", ok ? "PASS" : "FAIL");
	return ok ? 0 : 1;
}

static void usage(const char *p)
{
	fprintf(stderr,
		"usage: %s [--tmin T] [--tmax T] [--thmin TH] [--lo G] [--hi G]\n"
		"          [--n N] [--calls N] [--iters N] [--seed N] [--check]\n", p);
	exit(2);
}

int main(int argc, char **argv)
{
	double lo = 1, hi = 1000;
	long calls = 200000, seed = 0;
	int n = 7, iters = 10, check = 0, i;

	eemminput.mu = mpi;
	for (i = 1; i < argc; i++) {
		const char *a = argv[i];
		int has = i + 1 < argc;

		if (!strcmp(a, "--check")) check = 1;
		else if (!strcmp(a, "--tmin") && has) cut_tmin = atof(argv[++i]);
		else if (!strcmp(a, "--tmax") && has) cut_tmax = atof(argv[++i]);
		else if (!strcmp(a, "--thmin") && has) cut_thmin = atof(argv[++i]);
		else if (!strcmp(a, "--lo") && has) lo = atof(argv[++i]);
		else if (!strcmp(a, "--hi") && has) hi = atof(argv[++i]);
		else if (!strcmp(a, "--n") && has) n = atoi(argv[++i]);
		else if (!strcmp(a, "--calls") && has) calls = atol(argv[++i]);
		else if (!strcmp(a, "--iters") && has) iters = atoi(argv[++i]);
		else if (!strcmp(a, "--seed") && has) seed = atol(argv[++i]);
		else usage(argv[0]);
	}
	if (n < 1 || lo <= 0 || hi < lo || calls < 16 || iters < 1)
		usage(argv[0]);

	/* ranf(nonzero) seeds the generator on its FIRST call only
	   (kinc2/ranf.c); vegas always calls ranf(0) afterwards. */
	if (seed)
		(void)ranf(seed);

	if (check)
		return do_check(calls, iters);

	printf("# cuts: tmin=%g tmax=%g thmin=%g  calls=%ld iters=%d\n",
	       cut_tmin, cut_tmax, cut_thmin, calls, iters);
	printf("# sigma_raw is in raw/2s units, NOT pb (normalisation open)\n");
	printf("%10s %10s %14s %10s %6s %14s %12s\n", "sqrt_s", "Lambda",
	       "sigma_raw", "err", "chi2", "V_LL", "ratio");
	for (i = 0; i < n; i++) {
		double sq = grid_point(lo, hi, n, i), sig, err, chi2, v;

		run_point(sq, calls, iters, &sig, &err, &chi2);
		v = vll(sq * sq);
		printf("%10.4g %10.4f %14.6e %10.2e %6.2f %14.6e %12.6e\n",
		       sq, log(sq * sq / (mpi * mpi)), sig, err, chi2, v,
		       sig / v);
		fflush(stdout);
	}
	return 0;
}
