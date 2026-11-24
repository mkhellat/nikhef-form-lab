/*
 * Thin VEGAS integrator.
 *
 * Algorithm: G. P. Lepage, J. Comput. Phys. 27 (1978) 192.
 * This translation follows the course copy formcourse/kinc/vegas.c
 * (same text as 02S-sample-vegas-pipeline/vegas.c and the topic 03–05
 * copies). Histogram hooks are omitted. ranf / ipow / iipow are the
 * vendored kinc2 routines.
 *
 * The prototype matches kinc2/axolib.h so a driver can call vegas()
 * through that header without including this one.
 */
#ifndef VEGAS_H
#define VEGAS_H

void vegas(double (*fcn)(double *), double accuracy, int dimension,
	   long ncalls, int itmax, int printflag);

extern double avgi, sd, chi2a;

#endif
