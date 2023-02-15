#include <stdio.h>
#include <stdlib.h>

int iipow(int x,int y)
{
	int z, u;
/*
	determines x to the power y
*/
	if ( y == 0 ) return(1);
	if ( y < 0 ) {
		if ( x == 1 ) return(1);
		else if ( x == -1 ) {
			if ( ( y&1 ) != 0 ) return(-1);
			else return(1);
		}
		else if ( x == 0 ) {
			printf("Division by zero in iipow: %d^%d\n",x,y);
			exit(-1);
		}
		else return(0);
	}
	if ( y == 2 ) u = x*x;
	else {
		if ( ( y & 1 ) != 0 ) u = x;
		else u = 1;
		z = x;
		y >>= 1;
		while ( y ) {
			z = z*z;
			if ( ( y & 1 ) != 0 ) u *= z;
			y >>= 1;
		}
	}
	return(u);
}

