/*
	Main program for the reaction e-e+ -> mu-mu+

	Note that alfa(q^2) = alfa/(1-alfa/(3*pi)*log(-q^2/(A*m^2)))
		with A = exp(5/3)
*/
#include "axolib.h"
double mumass = (double)0.105658;
double emass  = (double)0.000511;
double alpha0 = ((double)(1.))/((double)(137.0368));
double pi;
/*
	Now set the energy of the beam in GeV:
*/
double Ebeam = 20.;

double alpha;
double s;
double overal;

double feemm(double *);

int main(int argc,char **argv)
{
	double sp;

	pi = acos(-1.);

	s  = 4. * Ebeam*Ebeam;
	// 4-momentum of the electron and positron in the Center of
	// Momentum frame would be p1=(Ebeam,pe) and p2=(Ebeam,-pe);
	//
	// 4-momentum oconservation and momemtum-energy relations:
	//    p1 + p2 = p3 + p4   (4 eqs)
	//    p1^2 = m1^2,
	//    p2^2 = m2^2,
	//    p3^2 = m3^2,
	//    p4^2 = m4^2         (4 eqs)
	//
	// There are a total of 16 4-momentum combinations (p_i +
	// p_j).(p_i + p_j) or equivalently p_i.p_j for any two-body
	// scatterings. Taking into account the above 8 equations, it
	// means only 8 of those combinations are independent.
	//
	// There are also 6 further implicit relations coming from the
	// commutative nature of inner product of 4-vectors:
	//
	//  p1.p2 = p2.p1, p1.p3 = p3.p1, ...
	//
	// This all would lead to ONLY 2 independant binary inner
	// products for the two-body scatterings or equivalently two
	// independent (p_i + p_j)^2 combinations.
	//
	// However, 3 mandelstam variables are conventionally used to
	// parameterize two-body scatterings.
	//
	// Mandelstam variables for two-body scattering would be:
	//
	//  t =  q^2 = (p1-p3)^2 = (p4-p2)^2  : 4-momentum-transfer
	//  s = -q^2 = (p1+p2)^2 = (p3+p4)^2  : total 4-momentum
	//  u        = (p1-p4)^2 = (p2-p3)^2
	//
	// s+t+u     = 2*me^2 + 2*mmu^2       : total mass present in
	//                                      the process
	//
	// In very high energies where we could basically ignore all
	// the masses present in the scattering, the Mandelstam
	// variables basically reduce to the inner product of all
	// 4-momenta:
	//
	//  t = -2*p1.p3 = -2*p2.p4
	//  s =  2*p1.p2 =  2*p3.p4
	//  u = -2*p1.p4 = -2*p2.p3
	//
       	alpha = alpha0/(1.-alpha0/(3.*pi)*log(s/(emass*emass*exp(5./3.))));
	//
	// QED effective coupling :
	//
	// at momentum transfers |q^2| > (30 GeV)^2 (which is our
	// case), one should take into account electroweak corrections
	//
	//   alpha(-q^2 = 0) (Thomson limit) ~ 1/137.036
	//   alpha(-q^2 = Mz^2)              ~ 1/128
	//
	// However, taking into account only the VP contribution from
	// electron loops, the effective qed coupling at one-loop
	// order for large momentum transfers (-q^2 >> m_e^2) would be
	// Peskin(7.96):
	//
	//  alpha_eff(-q^2) = alfa/(1-alfa/(3*pi)*log(-q^2/(A*m^2)))
	//  A = A = exp(5/3)
/*
	Next is the flux factor.
	Because the masses are identical we can take the root.
*/
	sp = 4. * ( Ebeam*Ebeam - emass*emass );
/*
		factors:
			19732.7^2 from 1/GeV^2 to picobarn
			1/(2*sp)  flux factor
			e^2 = 4*pi*alpha       twice
			1/4       spin averaging
*/
	overal = 19732.7*19732.7*(16.*pi*pi*alpha*alpha)/(2.*sp)/4.;
/*
	Now we initialize the histograms
	Note that 'instance' 1 is defined in the routine
	PlotBoundaries which can be found in the file boundaries.c
	Histograms will only be made if argc == 2. This means that the program
	is called with one argument. The argument should be the name of the
	output file.
*/
	if ( argc == 2 ) {
		if ( SetPlotFile(argc,argv) ) return(-1);
		if ( inplot(3) == 0 ) return(-1);
	}
	else if ( argc == 1 ) {
		printf("Correct syntax is: %s <nameoffilewithgraphs.tex>\n",*argv);
		exit(-1);
	}
/*
	The call to vegas. The arguments are:
		feemm   The function to be integrated
		0.0001  The required accuracy
		2       The number of dimensions
		1000    The number of points per iteration
		35      The maximum number of iterations
		1       Normal printout. Result of each iteration.
*/
	vegas(feemm,(double)0.00001,2,60000,6,1);
	return(0);
}
