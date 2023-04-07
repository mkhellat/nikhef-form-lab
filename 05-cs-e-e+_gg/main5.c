/*
	Main program for the reaction e-e+ -> gamma gamma
	Note that alfa(q^2) = alfa/(1-alfa/(3*pi)*log(-q^2/(A*m^2)))
		with A = exp(5/3)
*/
#include "axolib.h"
double emass = (double)0.000511;
double alpha0 = ((double)(1.))/((double)(137.0368));
double pi;

double Ebeam = 15.;

double anglecut = 0.05;

double alpha;
double s;
double overal;

double feegg(double *);

int main(int argc,char **argv)
{
	double sp;

	pi = acos(-1.);

	s  = 4. * Ebeam*Ebeam;

	alpha = alpha0/(1.-alpha0/(3.*pi)*log(s/(emass*emass*exp(5./3.))));

	sp = 4. * ( Ebeam*Ebeam - emass*emass );
/*
	Now we initialize the histograms
	Note that 'instance' 5 is defined in the routine
	PlotBoundaries which can be found in the file boundaries.c
	Histograms will only be made if argc == 2. This means that the program
	is called with one argument. The argument should be the name of the
	output file.
*/
	if ( argc == 2 ) {
		if ( SetPlotFile(argc,argv) ) return(-1);
		if ( inplot(5) == 0 ) return(-1);
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
*/
	overal = 19732.*19732.*(16.*pi*pi*alpha*alpha)/(2.*sp)/4.;

	vegas(feegg,(double)0.00000001,2,500000,10,1);
	return(0);
}
