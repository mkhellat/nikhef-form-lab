/*
 * ------------------------------------------------------------------
 *  mc-str-1d : Monte-Carlo integration of a 1-D function over the
 *              interval [0,1] using strata sampling
 * ------------------------------------------------------------------
 *  Test function fun1(x) is integrated over the Real interval 0 to 1.
 *  Random number generator is ranf(LONG dummy). 
 *
 *  Integration region is partitioned into np1 number of equisized
 *  intervals and in each interval np2 number of points are sampled
 *  with uniform distribution.
 *
 *  Sample total and its error is calculated for each stratum and are
 *  summed for the final result (the errors have to be added
 *  quadratically). To compare, the same integral is estimated without
 *  using stratified sampling.
 *
 */

#include "axolib.h"


double
fun1(double x)
{
  return(sinh(sin(x))/(pow(x,3) + 1));
}


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


int
main(int argc,char **argv)
{
  int n1, n2,
    n_strata = stoint(argv[1]),
    n_sample_points_per_stratum = stoint(argv[2]),
    n_total_sample_points = n_strata*n_sample_points_per_stratum;
  double x, y, fsum, ffsum, gsum, ggsum, hsum, hhsum, sigma2, error, error2, f;
  printf("number of strata: %d , number of points per stratum: %d\n", \
	 n_strata , n_sample_points_per_stratum);
  fsum = ffsum = 0;
  hsum = hhsum = 0;
  for ( n1 = 0; n1 < n_strata; n1++ ) {
    gsum = ggsum = 0;
    for ( n2 = 0; n2 < n_sample_points_per_stratum; n2++ ) {
      y = x = ranf(0);
      f = fun1(x);
      fsum = fsum + f;
      ffsum = ffsum + f*f;
      x = (y+n1)/n_strata;
      //printf("bin %d point no. %d --> %10.8f , f = %10.8f\n",n1+1,n2+1,x,f);
      f = fun1(x);
      gsum = gsum + f;
      ggsum = ggsum + f*f;
    }
    gsum = gsum / n_sample_points_per_stratum;
    ggsum = ggsum / n_sample_points_per_stratum;
    printf("...................... stratum %d ......................\n",n1+1);
    printf("  E(f) for stratum %d --> %10.8f\n",n1+1,gsum);
    printf("E(f^2) for stratum %d --> %10.8f\n",n1+1,ggsum);
    printf("E(f)^2 for stratum %d --> %10.8f\n",n1+1,gsum*gsum);
    sigma2 = ggsum - gsum*gsum;
    error2 = sigma2/(n_sample_points_per_stratum-1);
    printf("--- SAMPLE VARIANCE      for stratum %d : %10.8f\n",n1+1,sigma2);
    printf("--- SAMPLE STD DEVIATION for stratum %d : %10.8f\n",n1+1,sqrt(sigma2));
    printf("--- SAMPLE STD ERROR     for stratum %d : %10.8f\n",n1+1,sqrt(error2));
    /* (A) -----------------------------------------------------------
       Variance of the whole strata is the sum of the variances of all
       strata as each stratum is considered to be independent of any
       other stratum due to sample averages of all strata being
       independent of each other and leading to cancellation of
       terms.
         \bar{y_{st}} : average of all strata
	 \bar{y_h}    : average of stratum h
	 N_h          : number of sanple points in stratum h
	 N            : total number of sample points
	 \hat{T_{st}} : all strata total (sample)
	 \hat{T_h}    : stratum h total  (sample)

	          \hat{T_{st}} = \SIGMA_h \hat{T_h}         (A-1)

	 (A-1) => \bar{y_{st}} = \SIGMA_h N_h/N \bar{y_h}   (A-2)

	 (A-1) => VAR(\hat{T_{st}}) = VAR(\SIGMA_h \hat{T_h})
	          ;  after camnelation of terms
		  VAR(..) = \SIGMA_h VAR(\hat{T_h})         (A-3)
       
       (B) -----------------------------------------------------------       
       Now if the number of sample points per stratum are the same, it
       means that

          \frac{1}{n} VAR(str) =
	     \frac{1}{n} * \SIGMA_h {VAR(str_h)}            (B-1)
    */
    hhsum = hhsum + error2;
    /* (C) ----------------------------------------------------------- 
       to obtain the (estimated) total of the strata one should add
       the (estimated) total of each strata, i.e. (A-1)
       
          t_{str} = \SIGMA_h{t_h}.
	  
       And the estimated total for the sample points in the stratum h
       is t_h = N_h * \bar{y_h} = \frac{N_h}{n_h} \SIGMA_i{y_{hi}} and
       t_{str} = N * \bar{y_{str}} = \SIGMA_h{N_h * \bar{y_h}} which
       means \bar{y_{str}} = \SIGMA{ W_h * y_h} and W_h in our case is
       1/n_strata
    */
    hsum = hsum + gsum/n_strata;
  }

  printf("\n::::::::::::::::::::::: FINAL RESULT :::::::::::::::::::::::\n");
  fsum = fsum / n_total_sample_points;
  ffsum = ffsum / n_total_sample_points;
  sigma2 = ffsum - fsum*fsum;
  error = sqrt(sigma2/n_total_sample_points);
  printf(": For one range of %d points: %10.8f +/- %10.8f\n",n_total_sample_points,fsum,error);

  error = sqrt(hhsum)/n_strata;
  printf(": For %d ranges of each %d points: %10.8f +/- %10.8f\n",n_strata,n_sample_points_per_stratum,hsum,error);
  printf("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n\n");
  return(0);
}
