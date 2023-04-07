double ipow(double x,int y)
{
	int den;
	double z, u;
/*
	determines x to the power y
*/
	if ( y == 0 ) return(1.0);
	if ( y < 0 ) { den = 1; y = -y; }
	else den = 0;
	if ( y == 2 ) u = x*x;
	else {
		if ( ( y & 1 ) != 0 ) u = x;
		else u = 1.0;
		z = x;
		y >>= 1;
		while ( y ) {
			z = z*z;
			if ( ( y & 1 ) != 0 ) u *= z;
			y >>= 1;
		}
	}
	if ( den ) return(1.0/u);
	else return(u);
}

