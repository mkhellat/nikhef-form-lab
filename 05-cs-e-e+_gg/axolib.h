#ifndef AXOLIB__H
#define AXOLIB__H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct {
    double e;
    double px,py,pz;
    double m;
} PVEC;

#define DOTPR(p1,p2) (p1.e*p2.e-p1.px*p2.px-p1.py*p2.py-p1.pz*p2.pz)

typedef long LONG;
typedef double (*DFCN)(double *);

void vegas(
    DFCN fcn,          /* The function to be integrated */
    double accuracy,   /* The desired relative accuracy */
    int dimension,     /* Number of dimensions */
    long ncalls,       /* Number of points per iteration */
    int itmax,         /* Maximum number of iterations */
    int printflag);    /* How much output */

double ranf(LONG dummy);
double ipow(double,int);
int iipow(int,int);

double mgoto2(PVEC *, PVEC *, PVEC *, double *);
double mgoto3(PVEC *, PVEC *, PVEC *, PVEC *, double *);
 
extern double avgi,sd,chi2a;

/*
	The next variables are only needed when making graphs, also called 
	plots or histograms. They are single differential distributions.
	If no provisions are taken nothing will happen.
	To make histograms one should call inplot before calling vegas
	and fill the histograms inside the function that is called by vegas.
*/
 
#define LINEAR 1
#define LOGARITHMIC 2

typedef struct lplot {
	double *xvalue;
	double *average;
	int numxvalues;
	int numaverages;
} LPLOT;

typedef struct oneplot {
	char *title;
	char *latexname;
	double *xls;
	double *yls;
	int *nlsn;
	int *mlsn;
	double xmin;
	double xmax;
	double dls;
	double xlav;
	double xlsq;
	double xlava;
	double sxa;
	double top;
	double xltq;
	int nbins;
	int type;
} ONEPLOT;

typedef struct oneaverage {
	double zav;
	double yav;
	double zsv;
	double ysv;
	double ztv;
} ONEAVERAGE;

typedef struct lPlots {
	struct lPlots *next;
	ONEPLOT *theplots;
	ONEAVERAGE *theaverages;
	double fsqa;
	double vtot;
	int numlplots;
	int numaverages;
	int instance;
	int kt;
	int kk;
	int itt;
} LPLOTS;

typedef struct {/* Needed for communication between integration and inplot */
	double s1;
	double s2;
	double s3;
	double s4;
} NEEDEDINPLOT;

extern NEEDEDINPLOT neededinplot;
extern LPLOTS *currentplot;
extern LPLOT lplot;

LPLOTS *inplot(int instance);
void Replot(LPLOTS *lplots);
void xplot(LPLOTS *lplots,double value,double bin);
void PlotIt(LPLOTS *lplots,int now);
int PlotBoundaries(LPLOTS *lplots,int instance);
int SetPlotFile(int argc, char **argv);
int TestPlotBoundaries(LPLOTS *lplots,int instance);

#endif

