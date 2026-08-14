#include "axolib.h"

double rvec[100], seeds[24], carry = 0., twopm24;
static int mcall = 0;
static int ncall = 0;
void rcargo(LONG dummy);
void rcarry(int mcall);
static int notyet = 1;
static int i24 = 24, j24 = 10, iseeds[24];
static LONG icons = 2147483563;
static LONG ITWO24 = 0xFFFFFF;
static double TWOP12 = 4096.;

double
ranf(LONG dummy)
{
	if ( ncall == 0 ) {
		ncall = 1;
		if ( dummy ) rcargo(dummy);
	}
	if ( mcall == 0 ) {
		mcall = 100;
		rcarry(mcall);
	}
	return(rvec[--mcall]);
}

void rcargo(LONG inseed)
{
	LONG jseed = inseed, k;
	double twom24 = 1.0;
	int i;
	notyet = 0;
	for ( i = 0; i < 24; i++ ) {
		twom24 /= 2;
		k = jseed/53668;
		jseed = 40014*(jseed-k*53668) - k*12211;
		if ( jseed < 0 ) jseed += icons;
		iseeds[i] = jseed % ITWO24;
	}
	for ( i = 0; i < 24; i++ ) seeds[i] = iseeds[i]*twom24;
	i24 = 24;
	j24 = 10;
	carry = 0;
	if ( seeds[24-1] < seeds[14-1] ) carry = twom24;
	twopm24 = twom24;
}

void rcarry(int lenv)
{
	int ivec;
	double uni;
/*
	if not initialized -> default initialization with fixed number
*/
	if ( notyet ) rcargo((LONG)314159265);
	for ( ivec = 0; ivec < lenv; ivec++ ) {
		uni = seeds[--j24] - seeds[--i24] - carry;
		if ( uni < 0 ) { uni += 1.0; carry = twopm24; }
		else carry = 0.;
		seeds[i24] = uni;
		if ( i24 == 0 ) i24 = 24;
		if ( j24 == 0 ) j24 = 24;
		rvec[ivec] = uni;
	}
}

void rcarut(LONG *isdext)
{
	int i, icarry;
	for ( i = 0; i < 24; i++ ) isdext[i] = seeds[i]*TWOP12*TWOP12;
	icarry = ( carry > 0 );
	isdext[24] = 1000*j24 + 10*i24 + icarry;
}

void rcarin(LONG *isdext)
{
	int i;
	LONG isd;
	double twom24 = 1.0;
	for ( i = 0; i < 24; i++ ) twom24 /= 2;
	for ( i = 0; i < 24; i++ ) seeds[i] = isdext[i]*twom24;
	isd = isdext[24];
	carry = (isd % 10) * twom24;
	isd /= 10;
	i24 = isd % 100;
	j24 = isd / 100;
	twopm24 = twom24;
}
