#include "axolib.h"

/*
	Calculates the multiperipheral diagrams in the reactions
		ee->eemumu, ep->epmumu and pp->ppmumu.
	The parameters typeup and typedown determine which particles are involved.
	With the value 0 we have electrons.
	When protons are involved, they may be elastic (case 1)
	or inelastic (case 2).
	The up vertex is the one with p1 and p3 in the kinematics and the
	down vertex is the one with p2 and p5 in the vertex.
	The case that a proton is elastic uses the 'standard' elastic
	form factors of the proton. The case 2 gives the quasi-elastic formfactors
	which are a type of average in the resonance region. This requires also
	that we select a value (cq integrate over) the invariant mass of the
	outgoing hadron system. Hence the kinematics and the vegas call need
	one or two extra variables.
	Deep inelastic scattering is not treated here. The mass of the quarks
	would cause problems. There would also be some double counting with
	the quasi-elastic region. For deep inelastic scattering one should
	construct a separate program.

	The derivation of the matrix element is described in the paper.
	Once the Levi-Civita tensors have been computed properly in the kinematics
	routines (gamgam, orient, pickin) there should be no numerical problems.
*/

double peripp(int typeup, int typedown)
{
	double rho = 0.585;
	double u1,u2,xt,xu,tau,x,en,rhot,v1,v2,qqq,qdq,t11,t12,t21,t22,sum,ttt,x1;
	switch ( typeup ) {
		case 0:
			u1 = u2 = 1;
			break;
		case 1:
			xt = 1-extra.t1/.71;
			xt = xt*xt;
			xu = 2.79/xt;
			u1 = xu*xu;
			tau = extra.t1/(4*pickzz.w1);
			u2 = (1/(xt*xt)-xu*xu*tau)/(1-tau);
			break;
		case 2:
			x = extra.t1/(extra.t1-pickzz.w3);
			x1 = 1-x;
			x1 = x1*x1;
			x1 = x1*x1;
			en = pickzz.d1-extra.t1;
			tau = extra.t1/(4*pickzz.w1);
			rhot = rho-extra.t1;
			u1 = (-.86926*rho*rho*pickzz.d1/rhot/rhot
				-2.23422*pickzz.w1*x1/(x*(x*.96-1.26)+1))/extra.t1;
			u2 = (-tau*u1-.12549*pickzz.d1*extra.t1*rho/rhot/rhot*pickzz.d1*
				pickzz.d1/en/en/pickzz.w1)/(1-en*en/(4*pickzz.w1*extra.t1));
			break;
		default:
			printf("typeup = %d is not supported in peripp\n",typeup);
			exit(-1);
			break;
	}
	switch ( typedown ) {
		case 0:
			v1 = v2 = 1;
			break;
		case 1:
			xt = 1-extra.t2/.71;
			xt = xt*xt;
			xu = 2.79/xt;
			v1 = xu*xu;
			tau = extra.t2/(4*pickzz.w2);
			v2 = (1/(xt*xt)-xu*xu*tau)/(1-tau);
			break;
		case 2:
			x = extra.t2/(extra.t2-pickzz.w5);
			x1 = 1-x;
			x1 = x1*x1;
			x1 = x1*x1;
			en = pickzz.d2-extra.t2;
			tau = extra.t2/(4*pickzz.w2);
			rhot = rho-extra.t2;
			v1 = (-.86926*rho*rho*pickzz.d2/rhot/rhot
				-2.23422*pickzz.w2*x1/(x*(x*.96-1.26)+1))/extra.t2;
			v2 = (-tau*v1-.12549*pickzz.d2*extra.t2*rho/rhot/rhot*pickzz.d2*
				pickzz.d2/en/en/pickzz.w2)/(1-en*en/(4*pickzz.w2*extra.t2));
			break;
		default:
			printf("typedown = %d is not supported in peripp\n",typedown);
			exit(-1);
			break;
	}
	qqq = dotps.q1dq*dotps.q1dq;
	qdq = 4*dotps.w6-pickzz.w4;
	ttt = civita.epsi-levi.delta*(qdq+dotps.q1dq2);
	t22 = 512*(civita.bb*(levi.delta*levi.delta-levi.gram)-ttt*ttt
		-levi.sa1*civita.a6*civita.a6-levi.sa2*civita.a5*civita.a5
		-levi.sa1*levi.sa2*qqq);
	t12 = 128*(-civita.bb*(levi.dd2+civita.g6)
		-2*(extra.t1+2*dotps.w6)*(levi.sa2*qqq+civita.a6*civita.a6))*extra.t1;
	t21 = 128*(-civita.bb*(levi.dd4+civita.g5)
		-2*(extra.t2+2*dotps.w6)*(levi.sa1*qqq+civita.a5*civita.a5))*extra.t2;
	t11 = 64*(civita.bb*(qqq-levi.g4-qdq*(extra.t1+extra.t2+2*dotps.w6))
		-2*(extra.t1+2*dotps.w6)*(extra.t2+2*dotps.w6)*qqq)*extra.t1*extra.t2;
	sum = (((u1*v1*t11+u2*v1*t21+u1*v2*t12+u2*v2*t22)
		/(extra.t1*extra.t2*civita.bb))/(extra.t1*extra.t2*civita.bb))*.25;
	return(sum);
}


