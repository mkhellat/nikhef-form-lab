/*
 * Toy integrals for the thin VEGAS copy.
 *   ∫_0^1 4/(1+x^2) dx = π
 *   ∫_[0,1]^2 exp(x+y) dx dy = (e-1)^2
 * A perfectly flat integrand makes the iteration variance zero and
 * the 1/variance weight undefined, so the smoke stays curved.
 * Pass: PASS: vegas_toy
 */
#include "vegas.h"

#include <math.h>
#include <stdio.h>

static double f_pi(double *x)
{
	return 4.0 / (1.0 + x[0] * x[0]);
}

static double f_exp2(double *x)
{
	return exp(x[0] + x[1]);
}

static int check(const char *name, double exact, double rel)
{
	double err = fabs(avgi - exact) / fabs(exact);

	printf("# %s  avgi=%.8g  sd=%.3g  exact=%.8g  rel=%.3g\n", name, avgi,
	       sd, exact, err);
	if (!isfinite(avgi) || !isfinite(sd) || sd < 0.0 || err > rel) {
		fprintf(stderr, "FAIL: %s\n", name);
		return 1;
	}
	return 0;
}

int main(void)
{
	int bad = 0;

	vegas(f_pi, 1e-4, 1, 4000, 5, 0);
	bad |= check("pi", 3.141592653589793, 1e-3);

	vegas(f_exp2, 1e-4, 2, 4000, 5, 0);
	bad |= check("exp2", (exp(1.0) - 1.0) * (exp(1.0) - 1.0), 1e-3);

	if (bad)
		return 1;
	printf("PASS: vegas_toy\n");
	return 0;
}
