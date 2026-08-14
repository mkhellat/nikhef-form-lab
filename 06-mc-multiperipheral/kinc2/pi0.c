#include "axolib.h"

double pi0(int par)
{
	double sum, tt = extra.t1*extra.t2;
	double la = (eemminput.mu*eemminput.mu-extra.t1-extra.t2)*0.5;
	double part1 = 0,part2 = 0,part3 = 0,part4 = 0;
	la = -(la*la-tt);
	part1 = -64*levi.gram/(tt*tt);
	part2 = -16*extra.t1*levi.dd2/(tt*tt);
	part3 = -16*extra.t2*levi.dd4/(tt*tt);
	part4 = -4*tt*la/(tt*tt);
	sum = part1+part2+part3+part4;
	return(sum);
}


