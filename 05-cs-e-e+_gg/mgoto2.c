#include "axolib.h"

static int firstcall = 1;
static double pie;

double mgoto2(PVEC *p_in, PVEC *p1, PVEC *p2, double *y)
/*
    Routine for a two body decay. Note that there is no protection against
    special peaking behaviour near ct = +-1.
    Parameters:
        p_in: The incoming total 4-momentum
        p1:   The 4-momentum of the first outgoing particle.
        p2:   The 4-momentum of the second outgoing particle.
        y:    The array with the (two) random numbers between 0 and 1.
    The input is in p_in and y. Also the masses of p1 and p2 are input.
    The output is in p1 and p2.
    The return value of the routine is the value of the Jacobian of the
    transformation from the unit cube to the phase space.
*/
{
    double ch,sh,aa,a,pa,cp,sp,ct,st,e1p,p1a,ctp,stp,phi,p1x,p1y,p1z;

    if ( firstcall ) { firstcall = 0; pie = acos(-1); }
/*
    Determine the lorenz transformation that takes p_in from its rest frame
    to its actual values. Worry about special cases (division by zero....)
*/
    ch = p_in->e / p_in->m;
    aa = p_in->px*p_in->px+p_in->py*p_in->py;
    a  = sqrt(aa);
    pa = sqrt(aa+p_in->pz*p_in->pz);
    sh = pa / p_in->m;
    if ( pa == 0 ) {
        cp = ct = 1;
        sp = st = 0;
    }
    else {
        ct = p_in->pz/pa;
        if ( a == 0 ) {
            st = sp = 0; cp = 1;
        }
        else {
            st = a/pa;
            cp = p_in->px/a;
            sp = p_in->py/a;
        }
    }
/*
    In the CM-frame the vector p1 becomes
*/
    e1p = (p_in->m*p_in->m+p1->m*p1->m-p2->m*p2->m)/(2*p_in->m);
    p1a = sqrt((e1p-p1->m)*(e1p+p1->m));
    ctp = 1-2*y[0];
    stp = 2*sqrt(y[0]*(1-y[0]));
    phi = 2*pie*y[1];
    p1x = p1a*stp*cos(phi);
    p1y = p1a*stp*sin(phi);
    p1z = p1a*ctp;
/*
    Now we transform to the frame of p_in
*/
    p1->e  = ch*e1p+sh*p1z;
    p1->px = sh*cp*st*e1p+cp*ct*p1x-sp*p1y+cp*st*ch*p1z;
    p1->py = sh*sp*st*e1p+sp*ct*p1x+cp*p1y+sp*st*ch*p1z;
    p1->pz = sh*ct*e1p-st*p1x+ct*ch*p1z;
    p2->e = p_in->e-p1->e;
    p2->px = p_in->px-p1->px;
    p2->py = p_in->py-p1->py;
    p2->pz = p_in->pz-p1->pz;
/*
    Finally the jacobian of the transformation
*/
    return(p1a/(4*pie*p_in->m));
}
