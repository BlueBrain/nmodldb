VERBATIM
#include <math.h>
#include <time.h>
ENDVERBATIM

FUNCTION my_asin(arg) {
	VERBATIM
	double ret;
	ret=asin(*getarg(1));
	ENDVERBATIM
	my_asin=ret
}

FUNCTION my_sin(arg) {
	VERBATIM
	double ret;
	ret=sin(*getarg(1));
	ENDVERBATIM
	my_sin=ret
}

FUNCTION exp_i(arg) {
	VERBATIM
	double ret=1.0;
	double euler=exp(1.0);
	int n=0;
	for (n=0;n<(*getarg(1));n++) {
		ret*=euler;
	}
	ENDVERBATIM
	exp_i=ret
}

FUNCTION mytime() {
	VERBATIM
	double ret = 0.0;
	ret = (double)time(0)/3600.0;
	ENDVERBATIM
	mytime = ret
}

