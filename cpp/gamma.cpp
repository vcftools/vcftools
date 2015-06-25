#include "gamma.h"

double gammln(double xx)
{
	double x, y, tmp, ser;
	int j;
	static const double cof[14]={57.1562356658629235,-59.5979603554754912, 14.1360979747417471,-0.491913816097620199,.339946499848118887e-4, .465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3, -.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3, .844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5}; if (xx <= 0) throw("bad arg in gammln");
	y=x=xx;
	tmp = x+5.24218750000000000;
	tmp = (x+0.5)*log(tmp)-tmp;
	ser = 0.999999999999997092;
	for (j=0;j<14;j++) ser += cof[j]/++y;
	return tmp+log(2.5066282746310005*ser/x);
}


double gcf(double a, double x, double &gln)
{
	int i;
	double an,b,c,d,del,h;
	double EPS = numeric_limits<double>::epsilon();
	double FPMIN = numeric_limits<double>::min() / numeric_limits<double>::epsilon();
	gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;;i++) {
	an = -i*(i-a);
	b += 2.0;
	d=an*d+b;
	if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) <= EPS) break;
	}
	return exp(-x+a*log(x)-gln)*h;
}


double gser(double a, double x, double &gln)
{
	double sum,del,ap;
	gln=gammln(a);
	ap=a;
	del=sum=1.0/a;
	for (;;)
	{
		++ap;
		del *= x/ap;
		sum += del;
		if (fabs(del) < fabs(sum)*numeric_limits<double>::epsilon())
		{
			return sum*exp(-x+a*log(x)-gln);
		}
	}
	return 0;
}


double gammp(double a, double x)
{
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0 || (x != x) || (a != a))
	{
		return numeric_limits<double>::quiet_NaN();
	}
	if (x==0.0) return 0.0;
	if (x < (a+1.0))
	{
		gamser=gser(a,x,gln);
		return gamser;
	} else {
		gammcf = gcf(a,x,gln);
		return 1.0-gammcf;
	}
}

double gammq(double a, double x)
{
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0 || (x != x) || (a != a))
	{
		return numeric_limits<double>::quiet_NaN();
	}
	if (x == 0.0) return 1.0;
	if (x < (a+1.0)) {
		gamser=gser(a,x,gln);
		return 1.0-gamser;
	} else {
		gammcf = gcf(a,x,gln);
		return gammcf;
	}
}

	
