/*
 * ------------------------------------------------------------------
 *  mc-eqc-1d : Monte-Carlo integration of a 1-D function over the
 *              interval [0,1] using naive equalization of integral
 *              contribution of bins
 * ------------------------------------------------------------------
 * This program naively makes the contributions from different bins
 * close to each other by modifyig bin sizes and keeping the number of
 * sampled points per bin the same.
 *
 */

#include "axolib.h"

#define K 1000  // subintervals factor


//
// fun1: integrand
//
double
fun1(double x)
{
  return(3*pow(x,2));
}


//
// stoint: convert string to int for argumet parsing
//
int
stoint(char* str)
{
  int i = 0, val = 0, sign = 1;

  if ( str[0] == '-' ){
    sign = -1;
    i++;
  } 

  for (; str[i] != '\0'; i++ )
    val = val * 10 + str[i] - '0';

  return sign * val;
}


//
// main: takes number of bins and number of points per bin as
//       argument
//
int
main(int argc,char **argv)
{
  int n1, n2, itr = 0;
  // number of point per stratum should be set to 2; yet for
  // experimental purposes, user can play around with it.
  int n_bins = stoint(argv[1]),
    n_samples_per_bin = stoint(argv[2]), 
    n_total_samples = n_bins*n_samples_per_bin;
  int n_total_subintervals = K + n_bins,
    n_subintervals_per_bin = (int)((K / n_bins) + 1 + 0.5);
  double x, y, xsum, xran,
    fsum, ffsum,
    gsum, ggsum, hsum, hhsum, sigma2, error, error2, f;
  double Dx[n_bins]    // array of bin lengths
    , ic[n_bins];      // array of integral contribution for each bin
  int msum, m[n_bins]; // array of number of subintervals for each bin
  double dx[n_total_subintervals]; // array of subinterval lengths

  printf("number of bins: %d , number of points per bin: %d\n", \
	 n_bins , n_samples_per_bin);
  //
  // first array of uniforn DXs
  //
  for ( n1 = 0; n1 < n_bins; n1++ )
    Dx[n1] = 1.0 / n_bins;
  //
  // iterations
  //
  do{
    itr += 1;
    msum = 0;
    xran = 0;
    fsum = ffsum = 0;
    hsum = hhsum = 0;
    for ( n1 = 0; n1 < n_bins; n1++ ) {
      gsum = ggsum = 0;
      if ( n1 > 0 )
	xsum = xsum + Dx[n1-1];
      else
	xsum = 0;
      for ( n2 = 0; n2 < n_samples_per_bin; n2++ ) {
	y = x = ranf(0);
	f = fun1(x);
	fsum = fsum + f;
	ffsum = ffsum + f*f;
	xran = y * Dx[n1];
	x = xsum + xran;
	f = fun1(x);
	//printf("bin %d point no. %d --> %10.8f\n",n1+1,n2+1,x);
	gsum = gsum + f;
	ggsum = ggsum + f*f;
      }
      gsum = gsum / n_samples_per_bin;
      ggsum = ggsum / n_samples_per_bin;
      //
      // contribution of this bin to the integral and SEM
      //
      ic[n1] = gsum * Dx[n1];
      sigma2 = ggsum - gsum*gsum;
      error2 = sigma2/(n_samples_per_bin-1);
      //printf("--- bin %d DX: %f , I: %f , SEM squared : %f\n", n1+1,Dx[n1],ic[n1],sigma2); 
      hhsum = hhsum + error2;
      hsum = hsum + ic[n1];
    }
    error = sqrt(hhsum)/n_bins;
    //printf(": For %d ranges of each %d points: %10.8f +/- %10.8f\n",n_bins,n_samples_per_bin,hsum,error);
    printf(" iteration %d : %10.8f +/- %10.8f\n",itr,hsum,error);
    fsum = fsum / n_total_samples;
    ffsum = ffsum / n_total_samples;
    sigma2 = ffsum - fsum*fsum;
    error = sqrt(sigma2/n_total_samples);
    printf("   (no bins) : %10.8f +/- %10.8f\n",n_total_samples,fsum,error);
    //
    // evaluate new DXs
    //
    for (n1 = 0; n1 < n_bins; n1++){
      m[n1] = (int)((K * ic[n1] / hsum) + 1 + 0.5);
      //printf("--- m[%d] -> %d\n",n1,m[n1]);
      msum = msum + m[n1];
      for (n2 = msum - m[n1]; n2 < msum; n2++)
	dx[n2] = Dx[n1] / m[n1];
    }
    //printf(" --- total m -> %d\n", msum);
    for (n1 = 0; n1 < n_bins; n1++){
      Dx[n1] = 0;
      for (n2 = n1 * msum/n_bins; n2 < (n1 + 1) * msum/n_bins; n2++)
	Dx[n1] = Dx[n1] + dx[n2];
    }
    printf("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n");
  }
  while (itr < 100);

  return(0);
}
