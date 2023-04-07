/*
   #[ Definitions and variables :
*/

#define MAXLPLOTS 20
#define MAXAVERAGES 10
 
#include "axolib.h"

static LPLOTS *firstlplots = 0;

static char *nameoftheprogram = 0;
static char *nameofthefile = "histograms.tex";

LPLOTS *currentplot = 0;

/*
   #] Definitions and variables :
   #[ TestPlotBoundaries :
*/
 
int TestPlotBoundaries(LPLOTS *lplots,int instance)
{
	int error = 0;
	if ( lplots->numlplots > MAXLPLOTS ) {
		printf("Requested number of histograms for instance %d is greater than %d\n",instance,MAXLPLOTS);
		printf("Please readjust parameter MAXLPLOTS in file axolib.h\n");
		error = -1;
	}
	if ( lplots->numaverages > MAXAVERAGES ) {
		printf("Requested number of averages for instance %d is greater than %d\n",instance,MAXAVERAGES);
		printf("Please readjust parameter MAXAVERAGES in file axolib.h\n");
		error = -1;
	}
	return(error);
}

/*
   #] TestPlotBoundaries :
   #[ SetPlotFile :
*/

int SetPlotFile(int argc, char **argv)
{
	argc--;
	nameoftheprogram = *argv++;
	if ( argc != 1 ) {
		printf("Correct syntax is: %s <nameoffilewithgraphs>\n",nameoftheprogram);
		return(-1);
	}
	nameofthefile = *argv;
	return(0);
}

/*
   #] SetPlotFile :
   #[ DoInplotPrinting :

	In this routine we print the histograms of inplot.
	This is done in axodraw format to include the results in latex files.

*/

#define GRAPHWIDTH 300.0
#define GRAPHHEIGHT 400.0

static FILE *f = 0;

typedef struct {
	double top;
	double extra;
	int divisions;
	int subdivisions;
} TOPS;
static TOPS possibletops[] = {
			 {1.0,0.1,10,1}
			,{1.25,0.1,13,1}
			,{1.6,0.1,16,1}
			,{2.0,1.0,2,10}
			,{2.5,1.0,3,10}
			,{3.0,1.0,3,10}
			,{4.0,1.0,4,10}
			,{5.0,1.0,5,5}
			,{6.0,1.0,6,5}
			,{8.0,1.0,8,4}
			,{10.0,1.0,10,1}
		};

void DoInplotPrinting(LPLOTS *lplots)
{
	double yssq, maximum, av, max10, x, y, el1, el2, sxsq;
	int i, j, n, nn, m, k;
	ONEPLOT *opl;
	ONEAVERAGE *oav;
	if ( f == 0 ) {
		f = fopen(nameofthefile,"w");
		if ( f == 0 ) {
			printf("Cannot open output file %s for histograms\n",nameofthefile);
			exit(-1);
		}
	}
/*
	Write the header of the LaTeX file
*/
	fprintf(f,"\\documentclass[12pt]{article}\n");
	fprintf(f,"\\usepackage{a4}\n");
	fprintf(f,"\\usepackage{axodraw}\n");
	fprintf(f,"\\usepackage{graphicx}\n");
	fprintf(f,"\\usepackage{latexsym}\n");
	fprintf(f,"\\begin{document}\n");
	fprintf(f,"\\noindent\n");
/*
	Now one histogram at a time
*/
	for ( i = 0; i < lplots->numlplots; i++ ) {
		opl = lplots->theplots + i;
/*
		First determine the maximum value in the plot.
*/
		maximum = 0; av = 0;
		for ( j = 1; j <= opl->nbins; j++ ) {
			opl->xls[j] = opl->yls[j]/neededinplot.s1;
			if ( opl->xls[j] > maximum ) maximum = opl->xls[j];
			av = av + opl->xls[j];
		}
		av = av/opl->nbins;
/*
		Next we need the top of the figure such that the maximum is somewhere
		in the upper half of the plot, but if the figure is rather flat we
		want the average near the middle.
		And we like sensible units.
*/
		n = opl->nbins / 10;
		nn = opl->nbins % 10;
		if ( nn > 0 ) n++;
		if ( n <= 4 ) nn = 10;
		if ( n > 4 ) nn = 5;
		if ( n > 5 ) nn = 4;
		if ( n > 10 ) nn = 1;
		if ( maximum < av*2.1 ) maximum = av*2.1;
		if ( maximum != 0 ) {
			max10 = log10(maximum);
			if ( max10 >= 0 ) { m = max10; }
			else { m = max10-1.0; }
			max10 = max10 - m;
			max10 = pow(10,max10);
			for ( k = sizeof(possibletops)/sizeof(TOPS); k > 0; k-- ) {
				if ( max10 > possibletops[k-1].top ) { k--; break; }
			}
/*
			Now the real maximum is divisions*extra*pow(10,m)
*/
			maximum = possibletops[k+1].divisions*possibletops[k+1].extra
					*pow(10.0,m);
		}
		else {
			k = 0;
			maximum = 1.0;
		}

		fprintf(f,"\\begin{center}\n");
		fprintf(f,"\\begin{picture}(%6.2f,%6.2f)(0,0)\n"
					,GRAPHWIDTH+80,GRAPHHEIGHT+50);
		fprintf(f,"\\SetOffset(40,30)\n");
		fprintf(f,"\\SetColor{Black}\n");
		fprintf(f,"\\LinAxis(0,0)(%6.2f,0)(%d,%d,5,0,1.5)\n",GRAPHWIDTH,n,nn);
		fprintf(f,"\\LinAxis(%6.2f,%6.2f)(0,%6.2f)(%d,%d,5,0,1.5)\n"
				,GRAPHWIDTH,GRAPHHEIGHT,GRAPHHEIGHT,n,nn);
		fprintf(f,"\\LinAxis(0,0)(0,%6.2f)(%d,%d,-5,0,1.5)\n",GRAPHHEIGHT
				,possibletops[k].divisions,possibletops[k].subdivisions);
		fprintf(f,"\\LinAxis(%6.2f,0)(%6.2f,%6.2f)(%d,%d,5,0,1.5)\n"
				,GRAPHWIDTH,GRAPHWIDTH,GRAPHHEIGHT
				,possibletops[k].divisions,possibletops[k].subdivisions);
		fprintf(f,"\\rText(-25.0,%6.2f)[c][l]{$1/\\sigma\\ d\\sigma/d%s$}\n"
				,GRAPHHEIGHT/2.0,opl->latexname);
		for ( j = 1; j <= opl->nbins; j++ ) {
			x = (GRAPHWIDTH/opl->nbins)*(j-0.5);
			y = (opl->xls[j]/maximum)*GRAPHHEIGHT;
			fprintf(f,"\\Vertex(%6.2f,%6.2f){2.0}\n",x,y);
		}
		fprintf(f,"\\end{picture}\\\\\n\n");
		fprintf(f,"%s\\vspace{5mm}\\\\\n",opl->title);
/*
		Now the overflow and underflow bins
*/
		el1 = opl->yls[0]*opl->dls;
		el2 = el1/neededinplot.s1;
		fprintf(f,"Total underflow: %12.4e (= %12.4e \\%%) (%d points)\\\\\n",el1,100.*el2,opl->mlsn[0]);
		el1 = opl->yls[opl->nbins+1]*opl->dls;
		el2 = el1/neededinplot.s1;
		fprintf(f,"Total overflow:\\  %12.4e (= %12.4e \\%%) (%d points)\\\\\n",el1,100.*el2,opl->mlsn[opl->nbins+1]);
		sxsq = sqrt(opl->sxa/lplots->itt);
		fprintf(f,"Average value of this variable is %12.4e $\\pm$ %12.4e\\vspace{3mm}\\\\\n"
				,opl->xlava,sxsq);
		fprintf(f,"Total crosssection = %14.8g $\\pm$ %14.8g with $\\chi^2$ = %10.4g\\vspace{3mm}\\\\\n"
				,avgi,sd,chi2a);

		fprintf(f,"\\end{center}\n");		
		fprintf(f,"\\newpage\n");		
	}
/*
	Next the averages
*/
	if ( lplots->numaverages > 0 ) {
		fprintf(f,"The following are averages with error estimates\n");
		for ( i = 0; i < lplots->numaverages; i++ ) {
			oav = lplots->theaverages + i;
			yssq = sqrt(oav->ysv/lplots->itt);
			fprintf(f,"  %d:   %15.5e $\\pm$ %15.3e\n",i,oav->yav,yssq);
		}
	}
/*
	Finally the last part of the LaTeX file.
*/
	fprintf(f,"\\end{document}\n");
}

/*
   #] DoInplotPrinting :
   #[ inplot :
*/

LPLOTS *inplot(int instance)
{
	LPLOTS *lplots, *fl;
	ONEPLOT *opl;
	ONEAVERAGE *oav;
	int i, j, n, newinstance = 0;
/*
	First look whether we have this set of plots already
	If not, allocate a new one and put it in the linked list.
*/
	lplots = firstlplots; fl = 0;
	while ( lplots ) {
		if ( lplots->instance == instance ) goto inivars;
		fl = lplots;
		lplots = lplots->next;
	}
	newinstance = 1;
	lplots = (LPLOTS *)malloc(sizeof(LPLOTS));
	lplots->next = 0;
	lplots->numlplots = 0;
	lplots->instance = instance;
	if ( fl ) fl->next = lplots;
	else { firstlplots = lplots; }
/*
	Now we allocate the memory for the boundaries etc.
*/
	lplots->theplots = (ONEPLOT *)malloc(sizeof(ONEPLOT)*MAXLPLOTS);
	lplots->theaverages = (ONEAVERAGE *)malloc(sizeof(ONEAVERAGE)*MAXAVERAGES);
/*
	Now we have to read a routine with boundaries
*/
	if ( PlotBoundaries(lplots,instance) < 0 ) {
		return(0);
	}
/*
	Now initialize all other variables
*/
inivars:
	lplots->kt = 0;
	lplots->kk = 0;
	lplots->itt = 0;
	for ( i = 0; i < lplots->numlplots; i++ ) {
		opl = lplots->theplots + i;
		if ( opl->nbins < 1 ) opl->nbins = 1;
		n = opl->nbins + 2;
		if ( newinstance ) {
			opl->xls = (double *)malloc(sizeof(double)*2*n);
			opl->yls = opl->xls+n;
			opl->nlsn = (int *)malloc(sizeof(int)*2*n);
			opl->mlsn = opl->nlsn+n;
		}
		for ( j = 0; j < n; j++ ) {
			opl->yls[j] = 0.;
			opl->mlsn[j] = 0;
		}
		opl->dls = (opl->xmax-opl->xmin)/opl->nbins;
	}
	for ( i = 0; i < lplots->numaverages; i++ ) {
		oav = lplots->theaverages + i;
		oav->yav = oav->ysv = 0.;
	}
	if ( lplots->numlplots > lplot.numxvalues ) {
		if ( lplot.xvalue ) {
			free(lplot.xvalue);
			lplot.xvalue = 0;
		}
		lplot.xvalue = (double *)malloc(sizeof(double)*lplots->numlplots);
		lplot.numxvalues = lplots->numlplots;
	}
	if ( lplots->numaverages > lplot.numaverages ) {
		if ( lplot.average ) {
			free(lplot.average);
			lplot.average = 0;
		}
		lplot.average = (double *)malloc(sizeof(double)*lplots->numaverages);
		lplot.numaverages = lplots->numaverages;
	}
	currentplot = lplots;
	return(lplots);
}

/*
   #] inplot :
   #[ Replot :
*/

void Replot(LPLOTS *lplots)
{
	int i, j, n;
	ONEPLOT *opl;
	ONEAVERAGE *oav;
	lplots->kt++;
	lplots->fsqa = 0.;
	for ( i = 0; i < lplots->numlplots; i++ ) {
		opl = lplots->theplots + i;
		n = opl->nbins + 2;
		for ( j = 0; j < n; j++ ) {
			opl->xls[j] = 0.;
			opl->nlsn[j] = 0;
		}
		opl->xlav = 0.;
		opl->xltq = 0.;
		opl->xlsq = 0.;
	}
	for ( i = 0; i < lplots->numaverages; i++ ) {
		oav = lplots->theaverages + i;
		oav->zav = oav->zsv = 0.;
	}
}

/*
   #] Replot :
   #[ xplot :
*/

void xplot(LPLOTS *lplots,double value,double bin)
{
	int i, nlps;
	ONEPLOT *opl;
	ONEAVERAGE *oav;
	double aux;
	lplots->fsqa += value*value/bin;
	lplots->itt++;
	for ( i = 0; i < lplots->numlplots; i++ ) {
		opl = lplots->theplots + i;
		nlps = (lplot.xvalue[i]-opl->xmin)/opl->dls+1.;
		if ( nlps < 0 ) nlps = 0;
		if ( nlps > opl->nbins ) nlps = opl->nbins+1;
		opl->xls[nlps] += value/opl->dls;
		opl->nlsn[nlps]++;
		aux = value*lplot.xvalue[i];
		opl->xlav += aux;
		opl->xltq += value*aux/bin;
		opl->xlsq += aux*aux/bin;
	}
	for ( i = 0; i < lplots->numaverages; i++ ) {
		oav = lplots->theaverages + i;
		aux = value*lplot.average[i];
		oav->zav += aux;
		oav->ztv += value*aux/bin;
		oav->zsv += aux*aux/bin;
	}
}

/*
   #] xplot :
   #[ PlotIt :
*/

void PlotIt(LPLOTS *lplots,int now)
{
	double sxf, vbef, vu, al1, al2, sxt, sx2, xhelp;
	int i, j, nlps;
	ONEPLOT *opl;
	ONEAVERAGE *oav;
	if ( lplots->kk <= 0 ) {	/* First call */
		for ( i = 0; i < lplots->numlplots; i++ ) {
			opl = lplots->theplots + i;
			nlps = opl->nbins+2;
			for ( j = 0; j < nlps; j++ ) {
				opl->mlsn[j] = opl->nlsn[j];
				opl->yls[j] = opl->xls[j];
			}
		}
	}
	else {
		vbef = lplots->vtot;
		vu = neededinplot.s4/neededinplot.s3;
		vu = vu*vu;
		for ( i = 0; i < lplots->numlplots; i++ ) {
			opl = lplots->theplots + i;
			nlps = opl->nbins+2;
			for ( j = 0; j < nlps; j++ ) {
				if ( opl->nlsn[j] == 0 ) continue;
				if ( opl->mlsn[j] == 0 ) {
					opl->mlsn[j] = opl->nlsn[j];
					opl->yls[j] = opl->xls[j];
				}
				else {
					al1 = vu/opl->nlsn[j];
					al2 = vbef/opl->mlsn[j];
					opl->mlsn[j] = opl->mlsn[j]+opl->nlsn[j];
					opl->yls[j] = (al2*opl->xls[j]+al1*opl->yls[j])/(al1+al2);
				}
			}
		}
	}
	for ( i = 0; i < lplots->numlplots; i++ ) {
		opl = lplots->theplots + i;
		sxf = opl->xlsq - opl->xlav*opl->xlav;
		sxt = opl->xltq - opl->xlav*neededinplot.s3;
		sx2 = opl->xlsq/(opl->xlav*opl->xlav)
			 +lplots->fsqa/(neededinplot.s3*neededinplot.s3)
			 -2.*opl->xltq/(opl->xlav*neededinplot.s3);
		xhelp = opl->xlav/neededinplot.s3;
		sx2 = sx2*xhelp*xhelp;
		if ( lplots->kt == 1 ) {
	      opl->xlava = opl->xlav/neededinplot.s3;
    	  opl->sxa = sx2;
		}
		else if ( ( xhelp = sx2+opl->sxa ) != 0 ) {
			opl->xlava = (opl->xlav*opl->sxa/neededinplot.s3+opl->xlava*sx2)/xhelp;
			opl->sxa = opl->sxa*sx2/xhelp;
		}
	}
    lplots->vtot = neededinplot.s2/neededinplot.s1;
	lplots->vtot = lplots->vtot*lplots->vtot;
/*
	Now we do the averages
*/
	for ( i = 0; i < lplots->numaverages; i++ ) {
		oav = lplots->theaverages + i;
		sxf = oav->zsv - oav->zav*oav->zav;
		sxt = oav->zsv/(oav->zav*oav->zav)
			+lplots->fsqa/(neededinplot.s3*neededinplot.s3)
			-2.*oav->ztv/(oav->zav*neededinplot.s3);
		sx2 = oav->zav/neededinplot.s3;
		sx2 = sxt/(sx2*sx2);
		if ( lplots->kt != 1 ) {
			oav->yav = oav->zav/neededinplot.s3;
			oav->ysv = sx2;
		}
		else {
			if ( ( xhelp = sx2+oav->ysv ) != 0 ) {
				oav->yav = (oav->ysv*oav->zav/neededinplot.s3+oav->yav*sx2)/xhelp;
				oav->ysv = oav->ysv*sx2/xhelp;
			}
		}
	}
/*
	The actual printing
*/
	if ( now == 2 ) DoInplotPrinting(lplots);

	lplots->kk++;

}

/*
   #] PlotIt :
*/

