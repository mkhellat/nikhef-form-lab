/*
 * ------------------------------------------------------------------
 *  mc-vegas78-1d : Monte-Carlo integration of a 1-D function over the
 *                  interval [0,1] using a version of vegas algorithm
 * ------------------------------------------------------------------
 *
 * VEGAS sets the size of each bin in a way that the contributions of
 * sampled bins to the total sampled integral would get closer and
 * closer to each other.
 *
 */

#include "axolib.h"

#define K     1000  // subintervals factor for VEGAS
#define ITRS   100  // number of iterations


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
//       argument and perform vegas algorith until
//
int
main(int argc,char **argv)
{
  int n1, n2, itr = 0;
  int n_bins = stoint(argv[1]),
    n_samples_per_bin = stoint(argv[2]), 
    n_total_samples = n_bins*n_samples_per_bin;
  int n_total_sbincs = K + n_bins,
    n_sbincs_per_bin = (int)((K / n_bins) + 1 + 0.5);
  double x, y, xsum, xran,
    fsum, ffsum,
    gsum, ggsum, hsum, hhsum, sigma2, error, error2, f;
  double Dx[n_bins]    // array of bin lengths
    , ic[n_bins];      // array of integral contribution for each bin
  int msum, m[n_bins]; // array of number of sbincs for each bin
  double dx[n_total_sbincs]; // array of subinterval lengths
  double csigma[ITRS], I[ITRS], csigma_total, I_total; // comulative
						       // estimates

  printf("number of bins: %d , number of points per bin: %d\n", \
	 n_bins , n_samples_per_bin);
  //
  // first array of uniform DXs
  //
  for ( n1 = 0; n1 < n_bins; n1++ )
    Dx[n1] = 1.0 / n_bins;
  //
  // VEGAS iterations
  //
  do{
    itr += 1;
    msum = 0;
    xran = 0;
    fsum = ffsum = 0;
    hsum = hhsum = 0;
    error = 0;
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
      hhsum = hhsum + error2;
      hsum = hsum + ic[n1];
    }
    //
    // store this iteration results 
    //
    csigma[itr-1] = error = sqrt(hhsum) / n_bins;
    I[itr-1] = hsum;
    printf("iteration %d : %10.8f +/- %10.8f\n",itr,hsum,error);
    fsum = fsum / n_total_samples;
    ffsum = ffsum / n_total_samples;
    sigma2 = ffsum - fsum*fsum;
    error = sqrt(sigma2/n_total_samples);
    printf("cross check : %10.8f +/- %10.8f\n",n_total_samples,fsum,error);
    //
    // evaluate new DXs
    //
    for (n1 = 0; n1 < n_bins; n1++){
      m[n1] = (int)((K * ic[n1] / hsum) + 1 + 0.5);
      //
      // int alpha = 1; // values in the range [1,2] could be used 
      // double icr = ic[n1] / hsum;
      // m[n1] = (int)(K *pow((icr - 1)*(1 / log(icr)), alpha) + 1 + 0.5);
      // 
      // could be used instead for better performance and preventing
      // strange convergence behaviors.
      //      
      msum = msum + m[n1];
      //printf("     bin %d  m ,  ic , hsum : %d , %f , %f\n",n1+1,m[n1],ic[n1],hsum);
      for (n2 = msum - m[n1]; n2 < msum; n2++){
	dx[n2] = Dx[n1] / m[n1];
	//printf("--- dx[%d] : %f\n",n2,dx[n2]);
      }
    }
    printf(" --- total m -> %d\n", msum);
    for (n1 = 0; n1 < n_bins; n1++){
      Dx[n1] = 0;
      for (n2 = n1 * msum/n_bins; n2 < (n1 + 1) * msum/n_bins; n2++){
	Dx[n1] = Dx[n1] + dx[n2];
	//printf("+ %d + sum -> Dx[%d] = %f\n",n2,n1,Dx[n1]);
      }
      //printf("NEW Dx[%d] -> %f\n", n1, Dx[n1]);
    }
    printf("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n");
  }
  while (itr < ITRS);

  //
  // comulative result : for integrands with high narrow peaks,
  //                     comulative result should be estimated
  //                     differently
  //
  I_total = 0;
  csigma_total = 0;
  for ( n1 = 0; n1 < ITRS; n1++){
    I_total = I_total + (I[n1]/(csigma[n1]*csigma[n1]));
    csigma_total = csigma_total + (1/(csigma[n1]*csigma[n1]));   
  }
  csigma_total = sqrt(1/csigma_total);
  I_total = csigma_total * csigma_total * I_total;
  printf("************************************************************\n");
  printf("*\n");
  printf("* COMULATIVE RESULT : %10.8f +/- %10.8f\n", I_total, csigma_total);
  printf("*\n");
  printf("************************************************************\n");
  return(0);
}
