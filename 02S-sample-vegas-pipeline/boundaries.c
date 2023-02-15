/*
	This routine contains all definitions of sets of plot boundaries.
	If we keep the old ones this way we can just look them up. In the
	long run this may save much work.
*/
 
#include "axolib.h"

LPLOT lplot;
 
int PlotBoundaries(LPLOTS *lplots,int instance)
{
	static double pie = 0.;
	if ( pie == 0 ) pie = acos(-1.0);
	switch ( instance ) {
		case 0:
		case 1:
			break;
		case 2:
			lplots->numlplots = 2;
			lplots->numaverages = 0;
            if ( TestPlotBoundaries(lplots,instance) ) return(-1);
            lplots->theplots[0].xmin  =  0;
            lplots->theplots[0].xmax  =  1.;
            lplots->theplots[0].nbins = 40;
            lplots->theplots[0].title = strdup("x[0]");
            lplots->theplots[0].latexname = strdup("x_0");
            lplots->theplots[0].type = LINEAR;
            lplots->theplots[1].xmin  = 0.;
            lplots->theplots[1].xmax  = 1.0;
            lplots->theplots[1].nbins = 40;
            lplots->theplots[1].title = strdup("x[1]");
            lplots->theplots[1].latexname = strdup("x_1");
            lplots->theplots[1].type = LINEAR;
			break;
		case 3:
			lplots->numlplots = 2;
			lplots->numaverages = 0;
			if ( TestPlotBoundaries(lplots,instance) ) return(-1);
			lplots->theplots[0].xmin  = -1.;
			lplots->theplots[0].xmax  =  1.;
			lplots->theplots[0].nbins = 40;
			lplots->theplots[0].title = strdup("cos(theta-CM)");
			lplots->theplots[0].latexname = strdup("\\cos\\ (\\theta^{CM})");
			lplots->theplots[0].type = LINEAR;
			lplots->theplots[1].xmin  = 0.;
			lplots->theplots[1].xmax  = pie;
			lplots->theplots[1].nbins = 40;
			lplots->theplots[1].title = strdup("phi-CM");
			lplots->theplots[1].latexname = strdup("\\phi^{CM}");
			lplots->theplots[1].type = LINEAR;
			break;
		case 4:
			lplots->numlplots = 3;
			lplots->numaverages = 0;
			if ( TestPlotBoundaries(lplots,instance) ) return(-1);
			lplots->theplots[0].xmin  = -1.;
			lplots->theplots[0].xmax  =  1.;
			lplots->theplots[0].nbins = 40;
			lplots->theplots[0].title = strdup("cos(theta-CM)");
			lplots->theplots[0].latexname = strdup("\\cos\\ (\\theta^{CM})");
			lplots->theplots[0].type = LINEAR;
			lplots->theplots[1].xmin  = 0.;
			lplots->theplots[1].xmax  = pie;
			lplots->theplots[1].nbins = 40;
			lplots->theplots[1].title = strdup("phi-CM");
			lplots->theplots[1].latexname = strdup("\\phi^{CM}");
			lplots->theplots[1].type = LINEAR;
			lplots->theplots[2].xmin  =  0.;
			lplots->theplots[2].xmax  =  1.;
			lplots->theplots[2].nbins = 40;
			lplots->theplots[2].title = strdup("x[0]");
			lplots->theplots[2].latexname = strdup("x_0");
			lplots->theplots[2].type = LINEAR;
			break;
		case 5:
			break;
		case 6:
			lplots->numlplots = 4;
			lplots->numaverages = 0;
			if ( TestPlotBoundaries(lplots,instance) ) return(-1);
			lplots->theplots[0].xmin  = 0.;
			lplots->theplots[0].xmax  = 1.;
			lplots->theplots[0].nbins = 40;
			lplots->theplots[0].title = strdup("E-electron");
			lplots->theplots[0].latexname = strdup("E_e");
			lplots->theplots[0].type = LINEAR;
			lplots->theplots[1].xmin  = 0.;
			lplots->theplots[1].xmax  = 1.;
			lplots->theplots[1].nbins = 40;
			lplots->theplots[1].title = strdup("E-mu-plus");
			lplots->theplots[1].latexname = strdup("E_\\mu");
			lplots->theplots[1].type = LINEAR;
			lplots->theplots[2].xmin  = 0.;
			lplots->theplots[2].xmax  = 1.;
			lplots->theplots[2].nbins = 40;
			lplots->theplots[2].title = strdup("E-anti-nu-e");
			lplots->theplots[2].latexname = strdup("E_{\\overline{\\nu_e}}");
			lplots->theplots[2].type = LINEAR;
			lplots->theplots[3].xmin  = 0.;
			lplots->theplots[3].xmax  = 180.;
			lplots->theplots[3].nbins = 40;
			lplots->theplots[3].title = strdup("theta25");
			lplots->theplots[3].latexname = strdup("\\theta_{e\\mu}");
			lplots->theplots[3].type = LINEAR;
			break;
		default:
			printf("Error in PlotBoundaries: no boundaries found for plot instance %d\n",instance);
			return(-1);
			break;
	}
	return(0);
}
