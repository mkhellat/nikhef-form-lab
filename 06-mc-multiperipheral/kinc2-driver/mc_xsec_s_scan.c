/*
	Crude Monte Carlo scan of the (unnormalized-by-flux) integrated
	cross-section proxy integral( dj * |M|^2 ) d(unit cube) for
	e-e+->e-e+pi0 as a function of sqrt(s), to check directly whether the
	physical, phase-space-integrated cross section grows with s the way
	pi0.c's bare part1 does (see ../notes/topic5-s-scaling-pi0.md), or is
	tamed by the phase-space Jacobian/boundary the way the uam19 lecture
	promises ("only a behaviour like 1/(t1 t2) could give a decent
	crosssection").

	This is NOT a re-derivation of the flux-factor normalization from
	scratch -- that is read directly from the real reference driver,
	kinc2/eee.c (not vendored into this repo's kinc2/, since eee.c pulls
	in the VEGAS/plotting machinery this directory's tooling deliberately
	avoids -- see ../kinc2/README.md): eee.c's eemminput.overallconstant
	includes an explicit `/(2*s)` flux factor on top of couplings/hbarc^2
	that are s-independent. This driver reproduces that structure exactly
	(the `/(2*s)` division applied to the MC sample average below) without
	linking eee.c itself.

	Method: for each sqrt(s), draw N uniform points in the 4D unit cube
	(the same rannums orient()/pickin() consume), evaluate dj*pi0(0) at
	each accepted (kinematically allowed) point, and average -- this IS
	exactly what eepi()'s VEGAS integrand computes per call (see
	kinc2/eepi.c in the original archive), just without VEGAS's adaptive
	importance sampling, so treat the result as an order-of-magnitude
	Monte Carlo estimate, not a precision cross-section value.

	Result (see ../notes/topic5-s-scaling-pi0.md): the raw sample average
	of dj*pi0(0) itself grows roughly like s^1.2-1.5 (decreasing power as
	s increases, consistent with an s*log(s)-type growth rather than a
	hard power law) over the range sqrt(s) = 1..100 GeV. After dividing by
	the flux factor 2s from eee.c, the growth softens further but does
	NOT fully vanish at this sample size/precision -- i.e. even the full,
	phase-space-integrated proxy cross section is NOT flat in s the way a
	strict "decent 1/(t1 t2)" reading might suggest, at least not over an
	unrestricted (uncut, full 4pi) t1,t2 acceptance. This is left as an
	open numerical result rather than a fully closed derivation: real
	experiments always impose a minimum |t1|,|t2| acceptance cut (the
	electron/positron scattering angle cannot be exactly zero), which
	directly excludes the growing near-real-photon corner that part1's
	analytic s^2 growth (scripts/gram_s_scaling.py) lives in; deriving the
	exact asymptotic power/log law for the fully-inclusive, uncut
	integral is a genuine phase-space-integral calculation beyond what
	this Monte Carlo probe resolves.

	Depends on the vendored kinc2 sources in ../kinc2/ (see
	../kinc2/README.md), exactly like scan_pi0_scaling.c and
	check_gram_minors.c.

	Build (from this directory, kinc2-driver/):
		cc -Wall -g -I../kinc2 -c ../kinc2/orient.c ../kinc2/pickin.c \
			../kinc2/mapt.c ../kinc2/mapw.c ../kinc2/mapla.c \
			../kinc2/ipow.c ../kinc2/iipow.c ../kinc2/ranf.c \
			../kinc2/pi0.c
		cc -Wall -g -I../kinc2 -c mc_xsec_s_scan.c
		cc -o mc_xsec_s_scan *.o -lm
		./mc_xsec_s_scan
*/

#include "axolib.h"
#include <stdio.h>
#include <stdlib.h>

double pi0(int par);

EEMMINPUT eemminput;

static double me = 0.000511;   /* electron mass, GeV */
static double mpi = 0.135;     /* pi0 mass, GeV */

static double urand(void)
{
	return (double)rand() / ((double)RAND_MAX + 1.0);
}

int main(void)
{
	double sqrts_vals[] = {1.0, 3.0, 10.0, 30.0, 100.0};
	long N = 2000000;
	int i;
	long n;

	eemminput.mu = mpi;
	srand(12345);  /* fixed seed: reproducible across runs */

	printf("%10s %16s %16s %16s\n",
		"sqrt(s)", "raw_avg", "raw_avg/(2s)", "good/N");

	for (i = 0; i < (int)(sizeof(sqrts_vals) / sizeof(sqrts_vals[0])); i++) {
		double sqrts = sqrts_vals[i];
		double s = sqrts * sqrts;
		double sum = 0.0;
		long good = 0;

		for (n = 0; n < N; n++) {
			double rannums[4] = {urand(), urand(), urand(), urand()};
			double dj = orient(s, me, me, me, mpi, me, rannums, 0);
			double m2;

			if (dj == 0) continue;
			m2 = pi0(0);
			if (m2 < 0) continue;  /* eepi.c treats this as a bad point too */
			sum += dj * m2;
			good++;
		}

		printf("%10g %16.6e %16.6e %10ld/%-10ld\n",
			sqrts, sum / N, sum / N / (2 * s), good, N);
	}
	return 0;
}
