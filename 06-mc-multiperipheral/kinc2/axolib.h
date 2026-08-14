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

typedef struct {
    double re;  /* real part */
    double im;  /* imaginary part */
} COMPLEX;
typedef struct {
    double w1,w2,w3,w4,w5,d1,d2,d5,d7,sl1;
} PICKZZ;
typedef struct {
    double s1,s2,t1,t2;
} EXTRA;
typedef struct {
    double acc3, acc4;
} ACCURA;
typedef struct {
    double gram,dd1,dd2,dd3,dd4,dd5,delta,g4,sa1,sa2;
} LEVI;
typedef struct {
    double p12,p13,p14,p15,p23,p24,p25,p34,p35,p45,p1k2,p2k1;
} DOTP;
typedef struct {
    double e,e1,e2,e3,e4,e5,p,p3,p4,p5,ct3,st3,ct4,st4,ct5,st5,cp3,sp3,cp5,sp5;
} VARIAB;
typedef struct {
    double al3,al4,be4,be5,de3,de5,pp3,pp4,pp5;
} VARIAC;
typedef struct {
    double e6,e7,p6,p7,ct6,st6,ct7,st7,cp6,sp6,cp7,sp7,w,pp6,pp7;
} VARIAD;
typedef struct {
    double q1dq,q1dq2,w6,w7;
} DOTPS;
typedef struct {
    double epsi,g5,g6,a5,a6,bb;
} CIVITA;
typedef struct {
    double ctg,stg,cpg,spg;
} EXT;
typedef struct {
    double ctcm6,stcm6;
} ANGU;
typedef struct {
    double me, mu, mp, e1, e2, p1, p2, sq, overallconstant, ebeam;
} EEMMINPUT;
typedef struct {
    double me, mu, mp, e1, e2, p1, p2, sq, overallconstant, ebeam;
} PPMMINPUT;
 
extern PICKZZ pickzz;
extern EXTRA extra;
extern ACCURA accura;
extern LEVI levi;
extern DOTP dotp;
extern VARIAB variab;
extern VARIAC variac;
extern VARIAD variad;
extern CIVITA civita;
extern DOTPS dotps;
extern EXT ext;
extern ANGU angu;
extern EEMMINPUT eemminput;
extern double qve[4];

double mapt(double,double,double *,double,int);
double mapw(double,double,double *,double,int);
double mapla(double smin, double smax, double y, double z
                    , double *ds, double x);
double pickin(double s, double m1, double m2, double m3, double m4, double m5
            ,double *rannums, int option);
double orient(double s, double m1, double m2, double m3, double m4, double m5
            ,double *rannums, int options);
double gamgam(double ebeam, double m1, double m2, double m3, double m5,
            double m6, double m7, double m4min, double m4max,
            double *rannums, int options, int nm);
double peripp(int typeup, int typedown);
double pi0(int par);

double ppmm(double *);
double epmm(double *);
double eemm(double *);
double eepi(double *);

#endif

