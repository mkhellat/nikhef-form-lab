#include "axolib.h"

extern double s;
extern double pi;
extern double emass;
extern double mumass;
extern double overal;
extern double Ebeam;

double meemm(double);

double feemm(double *x)
/*
    Two dimensional kinematics function for e-e+ -> mu-mu+
    The argument x contains the (two) random numbers between 0 and 1.
    The return value is the product of the matrix element squared,
    all kinematical Jacobians and the overal constants. These are all
    relevant factors with the exception of the weights used in the
    integration program.
*/
{
    PVEC p_in, pa, p1, p2;
    double jacobian, t, total;
/*
        First compose the incoming state. We use p_in, the sum of the incoming
        4-vectors and pa, the 4-vector of the electron. We don't need pb.
*/
    p_in.e = 2*Ebeam;      pa.e = Ebeam;
    p_in.px = 0.;          pa.px = 0;
    p_in.py = 0.;          pa.py = 0;
    p_in.pz = 0.;          pa.pz = sqrt((Ebeam-emass)*(Ebeam+emass));
    p_in.m = sqrt(s);      pa.m = emass;
/*
        Next we need to provide the masses of the outgoing particles.
        They are needed by mgoto2.
*/
    p1.m = p2.m = mumass ;
/*
        Now the call to mgoto2. In principle it should return a nonzero
        value unless we are exactly at the edge of phase space.
*/
    if ( ( jacobian = mgoto2(&p_in,&p1,&p2,x) ) == 0 ) return(0);
/*
        Convert the values of the 4-momenta into variables that we like more
*/
    t = mumass*mumass+emass*emass-2.*DOTPR(pa,p1);
/*
        Put now cuts if needed. If we reject a point because it is outside
        the acceptance of the detector, there is no need to evaluate the
        matrix element.

    if ( outside_acceptance ) return(0);
*/
    //
    // |cos(\theta_{CM})| >= 0.8 not accepted
    //
    if ( (1-2*x[0]) >= 0.8 || (1-2*x[0]) <= -0.8 ) return(0);
    //
    // DeadZone : 10 < \theta_{CM} < 25 AND 17 < \phi_{CM} < 39
    //
    if ( (1-2*x[0]) < cos(25) && (1-2*x[0]) > cos(10) && (2*pi*x[1]) < (39*pi/180) && (2*pi*x[1]) > (17*pi/180) ) return(0);
/*
    Next we fill the histograms.
    We have defined 2 histograms in PlotBoundaries.
    The boundaries are given in 'instance' = 3.
*/
    lplot.xvalue[0] = 1-2*x[0];
    lplot.xvalue[1] = 2*pi*x[1];
/*
        Evaluate the matrix element and compose the value of the output
*/
    total = overal*jacobian*meemm(t);
    return(total);
}
