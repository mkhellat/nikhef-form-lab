#include "axolib.h"

#define MAXDIM 20
#define MAXDIV 100

int ndo,it;
double si,si2,swgt,schi,scalls;
double xi[MAXDIV+1][MAXDIM],d[MAXDIV+1][MAXDIM],di[MAXDIV+1][MAXDIM];
int nxi[MAXDIV+1][MAXDIM];
double avgi,sd,chi2a;

NEEDEDINPLOT neededinplot;

void
vegas(
	DFCN fcn,          /* The function to be integrated */
	double accuracy,   /* The desired relative accuracy */
	int dimension,     /* Number of dimensions */
	long ncalls,       /* Number of points per iteration */
	int itmax,         /* Maximum number of iterations */
	int printflag)     /* How much output */
{
	int i,j,k,ndmx,nd,ng,npg,ndm,iaj,iaj1,mds;
	double calls,dxg,dv2g,xnd,xjac,rc,dr,xn,xo,ti,tsi,fb,f2b,wgt,f,f2,ti2;
	double xin[MAXDIV+1],r[MAXDIV+1],dx[MAXDIM],dt[MAXDIM];
	int ia[MAXDIM],kg[MAXDIM];
	double xl[MAXDIM],xu[MAXDIM],qran[MAXDIM],x[MAXDIM];
	double alpha = 1.5, one = 1.0;
	int now = 0;

	ndmx = MAXDIV;
	if ( ndmx*ndmx > ncalls ) ndmx = sqrt(((double)1.0)*ncalls);

	mds = 1;
	for ( i = 0; i < dimension; i++ ) { xl[i] = 0.; xu[i] = one; }
	ndo = 1;
	for ( j = 0; j < dimension; j++ ) xi[0][j] = one;

	it = 0;
	scalls = schi = swgt = si2 = si = 0;

	nd = ndmx;
	ng = 1;
	if ( mds != 0 ) {
	    ng = pow( ncalls*0.5, 1.0/dimension );
    	mds = 1;
		if ( (2*ng-ndmx) >= 0 ) {
			mds = -1;
			npg = ng/ndmx+1;
			nd = ng/npg;
			ng = npg*nd;
		}
	}
	k = iipow(ng,dimension);
	npg = ncalls/k;
	if ( npg < 2 ) npg = 2;
	calls = npg*k;
	dxg = one/ng;
	dv2g = ipow(dxg,2*dimension)/(npg*npg*(npg-one));
	xnd = nd;
	ndm = nd-1;
	dxg *= xnd;
	xjac = one;
	for ( j = 0; j < dimension; j++ ) {
		dx[j] = xu[j]-xl[j];
		xjac *= dx[j];
	}
	if ( nd != ndo ) {
		rc = ndo/xnd;
		for ( j = 0; j < dimension; j++ ) {
			i = k = 0;
			dr = xn = xo = 0;
			while ( i < ndm ) {
				while ( rc > dr ) {
					dr += one;
					xo = xn;
					xn = xi[k++][j];
				}
				dr -= rc;
				xin[i++] = xn-(xn-xo)*dr;
			}
			for ( i = 0; i < ndm; i++ ) xi[i][j] = xin[i];
			xi[nd-1][j] = one;
		}
		ndo = nd;
	}
	if ( printflag ) {
		if ( printflag == 10 ) {
			printf("vegas: dimension = %d, ncalls = %8.0f, itmax = %d",dimension,calls,itmax);
			printf(", acc = %e, mds = %d, nd = %d\n",accuracy,mds,nd);
		}
		else {
			printf("input parameters for vegas:  dimension = %3d  ncalls = %8.0f\n",dimension,calls);
			printf("                             it = %5d  itmax =    %5d\n",it,itmax);
			printf("                             acc = %e\n",accuracy);
			printf("                             mds = %3d\n",mds);
			printf("                             nd = %4d\n",nd);
		}
	}
reiterate:
	it++;
	tsi = ti = 0.;
	if ( currentplot ) {
		now = 0;
		Replot(currentplot);
	}
	for ( j = 0; j < dimension; j++ ) {
		kg[j] = 1;
		for ( i = 0; i < nd; i++ ) {
			nxi[i][j] = 0;
			d[i][j] = di[i][j] = ti;
		}
	}
reloop:
	f2b = fb = 0;
	for ( k = 0; k < npg; k++ ) {
		for ( j = 0; j < dimension; j++ ) qran[j] = ranf(0);
		wgt = xjac;
		for ( j = 0; j < dimension; j++ ) {
			xn = (kg[j]-qran[j])*dxg;
			ia[j] = xn;
			iaj   = ia[j];
			iaj1  = iaj-1;
			if ( iaj > 0 ) {
				xo = xi[iaj][j]-xi[iaj1][j];
				rc = xi[iaj1][j]+(xn-iaj)*xo;
			}
			else {
				xo = xi[iaj][j];
				rc = (xn-iaj)*xo;
			}
			x[j] = xl[j]+rc*dx[j];
			wgt *= xo*xnd;
		}
		f    = fcn(x)*wgt;
		if ( currentplot ) {
			xplot(currentplot,f/calls,wgt/calls);
		}
		f2   = f*f;
		fb   += f;
		f2b  += f2;
		for ( j = 0; j < dimension; j++ ) {
			iaj  = ia[j];
			(nxi[iaj][j])++;
			di[iaj][j] += f/calls;
			if ( mds >= 0 ) d[iaj][j] += f2;
		}
	}
	f2b = f2b*npg;
	f2b -= fb*fb;
	ti  += fb;
	tsi += f2b;
	if ( mds < 0 ) {
		for ( j = 0; j < dimension; j++ ) d[ia[j]][j] += f2b;
	}
	for ( k = dimension-1; k >= 0; k-- ) {
		kg[k] = ( kg[k] % ng ) + 1;
		if ( kg[k] != 1 ) goto reloop;
	}
	ti /= calls;
    tsi *= dv2g;
	ti2 = ti*ti;
	wgt = ti2/tsi;
	si += ti*wgt;
	si2 += ti2;
	swgt += wgt;
	schi += ti2*wgt;
	scalls += calls;
	avgi = si/swgt;
	sd = swgt*it/si2;
	chi2a = 0;
	if ( it > 1 ) chi2a = sd*(schi/swgt-avgi*avgi)/(it-1);
	sd = sqrt(one/sd);
	if ( printflag ) {
		tsi = sqrt(tsi);
		if ( printflag == 10 ) {
			printf("%d%20.8g%12.4g%20.8g%12.4g%12.4g\n",it,ti,tsi,avgi,sd,chi2a);
		}
		else {
			printf("\nintegration by vegas.\n\n");
			printf("iteration no %3d.      integral =%14.8g\n",it,ti);
			printf("                       std dev  =%14.8g\n",tsi);
			printf("accumulated results.   integral =%14.8g\n",avgi);
			printf("                       std dev  =%14.8g\n",sd);
			printf("                 chi^2 per itn  =%10.4g\n",chi2a);
		}
		if ( printflag < 0 ) {
			for ( j = 0; j < dimension; j++ ) {
				printf("data for axis%2d\n",j+1);
				printf("       x         delt i     convce        ");
/*				printf("       x         delt i     convce        "); */
				printf("       x         delt i     convce\n");
				for ( i = 0; i < nd; i++ ) {
					printf("%12.4g%12.4g%12.4g",xi[i][j],di[i][j],d[i][j]);
/*					if ( (i%3) == 2 ) printf("\n"); */
					if ( (i%2) == 1 ) printf("\n");
					else              printf("     ");
				}
/*				if ( (nd%3) != 0 ) printf("\n"); */
				if ( (nd%2) != 0 ) printf("\n");
			}
		}
	}
	if ( currentplot ) {
		if ( ( fabs(sd/avgi) <= fabs(accuracy) ) || ( it >= itmax ) ) now = 2;
		neededinplot.s1 = avgi;
		neededinplot.s2 = sd;
		neededinplot.s3 = ti;
		neededinplot.s4 = tsi;
        PlotIt(currentplot,now);
	}
	for ( j = 0; j < dimension; j++ ) {
		dt[j] = 0;
		for ( i = 0; i < nd; i++ ) {
			if ( nxi[i][j] > 0 ) d[i][j] /= nxi[i][j];
			dt[j] += d[i][j];
		}
	}
	for ( j = 0; j < dimension; j++ ) {
		rc = 0;
		for ( i = 0; i < nd; i++ ) {
			r[i] = 0;
			if ( d[i][j] > 0 ) {
				xo = dt[j]/d[i][j];
				r[i] = pow((xo-one)/(xo*log(xo)),alpha);
			}
			rc += r[i];
		}
		rc /= xnd;
		k = 0;
		dr = r[k];
		xo = 0;
		xn = xi[k][j];
		for ( i = 0; i < ndm; ) {
			while ( rc > dr ) {
				k++;
				dr += r[k];
				xo = xn;
				xn = xi[k][j];
			}
			dr -= rc;
			xin[i++] = xn-(xn-xo)*dr/r[k];
		}
		for ( i = 0; i < ndm; i++ ) xi[i][j] = xin[i];
		xi[nd-1][j] = one;
	}
	if ( ( it < itmax ) && ( fabs(accuracy) < fabs(sd/avgi) ) ) goto reiterate;
}
