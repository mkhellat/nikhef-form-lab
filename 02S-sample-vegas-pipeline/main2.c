/*
	Main program for the integration of the test function fun2
	It is a non-trivial two dimensional function
	When we activate the object WITHCUT a rather nontrivial piece
	it cut out from the integration region.

*/
//#define WITHCUT

#include "axolib.h"

int dimension = 2;

double fun2(double *x)
{
  double total;
#ifdef WITHCUT
  double y;
#endif
  total = 1.0/ipow(x[0]+x[1]*x[1]+0.1+ x[0]*x[1],2);
  lplot.xvalue[0] = x[0];
  lplot.xvalue[1] = x[1];
#ifdef WITHCUT
  y = x[0]*(1-x[0])+x[0]*x[1]*0.33+0.5*(x[1]-x[0]);
  if ( y > 0.2 && y < 0.3 ) total = 0;
#endif
  return(total);
}

int main(int argc,char **argv)
{
/*
    We initialize the histograms
    Note that 'instance' 1 is defined in the routine
    PlotBoundaries which can be found in the file boundaries.c
    Histograms will only be made if argc == 2. This means that the program
    is called with one argument. The argument should be the name of the
    output file.

    Here is what happens:
    
    (1) SetPlotFile() is called -> declared in "axolib.h" and defined in "inplot.c"
        (1-a) the fucntion fills in the value for __nameofthefile__
              and __nameoftheprogram__
	
    (2) vegas() calls PlotIT() and RePlot()
        (2-a) DoInplotPrinting() is called by the very two functions
              which would fill in the tex code using fprintf() to the
              filedescriptor for __nameofthefile__

   FURTHER NOTES: iipow() and ipow() functions are used by vegas and
   in the definition of fun2() whose definitions reside inside ipow.c
   amd iipow.c
*/
  {
    if ( SetPlotFile(argc,argv) ) return(-1);
    if ( inplot(2) == 0 ) return(-1);
  }

  vegas(fun2,(double)0.000001,dimension,8000,10,1);
  return(0);
}
