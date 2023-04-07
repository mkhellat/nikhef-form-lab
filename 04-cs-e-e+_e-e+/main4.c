/*
	Main program for the reaction e-e+ -> e-e+
	Note that alfa(q^2) = alfa/(1-alfa/(3*pi)*log(-q^2/(A*m^2)))
		with A = exp(5/3)
*/
#include "axolib.h"
double emass = (double)0.000511;
double alpha0 = ((double)(1.))/((double)(137.0368));
double pi;

double Ebeam = 15.;

// this is the anglecut var that is to be used in the kinematic file
// feeee.c to cut out the infinite part of the t-channel graph where
// clearly no detector is going to be present.
// 
// In the Center of Momentum frame:
//     t = - (s − 4*me*me)/2 * (1 - cos(theta))
//
// Matrix element expression from  meeee.c is:
//     48 + 16/s2*t2 + 32/s*t + 32*s/t + 16*s2/t2 - 64*me2/s2*t
//     - 64*me2*s/t2 + 64*me4/s2 - 64*me4/s/t + 64*me4/t2
//
//   p_a , p_b -> p_1 , p_2
//   a: incident beam
//   b: target particle
//   1 , 2: outgoing particles
//   d\sigma = P_{fi} *
//            ( V^2 / 2 * E_a * E_b * |v| ) *
//	      [ ( V / 2 * \pi^3 ) * ( d^3 p_1 / 2 * E_1 ) ] *
//	      [ ( V / 2 * \pi^3 ) * ( d^3 p_2 / 2 * E_2 ) ]
//
//            ::: transition rate per unit volume :::
//            P_{fi} = | A_{fi} | ^ 2 / V * T
//
//            ::: matrix element :::
//            A_{fi} = -i * \int {d^4 x} \Psi_3(x) \hat{V} \Psi_1(x)
//	             = -i * \int {d^4 x} j_{\mu} * A^{\mu}
//
//            ::: arbitrary normalization :::
//            \rho_i = 2 * E_i / V
//
//            ::: incident flux factor :::
//	      |v| * 2 * E_a / V
//
//            ::: number of target particles per unit volume :::
//	      2 * E_b / V
//
//            ::: final two-body phase space factor for integration :::
//	      [ ( V / 2 * \pi^3 ) * ( d^3 p_1 / 2 * E_1 ) ] *
//	      [ ( V / 2 * \pi^3 ) * ( d^3 p_2 / 2 * E_2 ) ]

double anglecut = 0.05;

double alpha;
double s;
double overal;

double feeee(double *);

int main(int argc,char **argv)
{
	double sp;

	pi = acos(-1.);

	s  = 4. * Ebeam*Ebeam;

	alpha = alpha0/(1.-alpha0/(3.*pi)*log(s/(emass*emass*exp(5./3.))));

	sp = 4. * ( Ebeam*Ebeam - emass*emass );
/*
	Now we initialize the histograms
	Note that 'instance' 2 is defined in the routine
	PlotBoundaries which can be found in the file boundaries.c
	Histograms will only be made if argc == 2. This means that the program
	is called with one argument. The argument should be the name of the
	output file.
*/
	if ( argc == 2 ) {
		if ( SetPlotFile(argc,argv) ) return(-1);
		if ( inplot(4) == 0 ) return(-1);
	}
	else if ( argc == 1 ) {
		printf("Correct syntax is: %s <nameoffilewithgraphs.tex>\n",*argv);
		exit(-1);
	}
/*
		factors:
			e^2 = 4*pi*alpha       twice
			1/(2*sp)  flux factor
			1/4       spin averaging
			19732     from GeV to picobarn

	::: unit conversion :::
	1 barn    10^-28 m^2, 100 fm^2
	1 pb      10^-40 m^2
	1 \bar{h}^2 c^2 / GeV^2 == 0.3894         mb
	                        == 3.8935 * 10^-4  b
				== 3.8935 * 10^+8 pb
				== (1.9732)^2 * (10^4)^2

        ::: flux factor :::
	In the CM frame
	
	    pa = (Ebeam , p) ; pb = (Ebeam, -p)
	    
        which means that incident beam momentum is 2 * p
	
	    Ebeam * Ebeam - p * p = emass * emass
	    
	In other words, the incident beam momentum squared is
	
	    4 * p * p = 4 * (Ebeam * Ebeam - emass * emass)
	              = sp

        Flux of the incident beam is the number of incident particles
        passing through unit area per unit time; lets assume we have
        $N_a$ particles of incident beam particles per unit volume and
        $N_b$ target particles per unit volume. This means that the
        incident beam flux is simply $|v| * N_a$ and our cross section
        should be per incident flux per target particle number
        density, i.e. we should divide matrix element by

	     |v| * N_a * N_b

        Taking into account that the matrix element itself contains
        N_a * N_b * N_1 * N_2, finally what we call __flux factor__
        would be

	      4 * E_a * E_b * |v| =
	         2 *
		 sqrt (
		    [ s - (m_a + m_b) * (m_a + m_b) ] *
		    [ s - (m_a - m_b) * (m_a - m_b) ]
		      )
		 = 2 * sqrt( sp * s ) ~ 2 * sp

*/
	overal = 19732.*19732.*(16.*pi*pi*alpha*alpha)/(2.*sp)/4.;

	vegas(feeee,(double)0.000001,2,500000,10,1);
	return(0);
}
