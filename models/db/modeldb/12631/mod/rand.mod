: $Id: rand.mod,v 1.4 1995/08/23 18:33:27 billl Exp $
COMMENT
$Header: /usr/site/nrniv/local/mod/RCS/rand.mod,v 1.4 1995/08/23 18:33:27 billl Exp $

Author: Stephen Fisher
Date:   December 1992
Email:  fisher@james.psych.yale.edu

Misc. random routines:
	fseed(seed)
		- set seed
		- return seed

	n_rand()
		- uniform distribution (0.0 <= rand < 1.0)

	fran(low, high)
		- returns random number between low and high

	u_rand()
		- uniform distribution (0.0 <= rand <= 1.0)

	norm()
		- gaussian distribution around 0

	pois(mean)
		- poisson distribution

	poisrand(mean)
		- integer poisson distribution using scop

** Note that with a SUFFIX equal to "nothing" these functions do not
have a suffix in hoc.  Thus to call norm() in hoc use simply type
"norm()" <- without the quotes.
ENDCOMMENT
					       
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX nothing
}

VERBATIM
#include <stdlib.h>
#include <math.h>
/* #include <values.h> /* contains MAXLONG */
#include <limits.h> /* contains MAXLONG */
/* Michael Hines fix for cygwin on mswin */
#if !defined(MAXLONG)

#define MAXLONG LONG_MAX
#endif
/*
extern double drand48();
#define random()                drand48()*MAXLONG
#define initstate(c1,c2,c3)     srand48(c1)
*/

static long state2[32] = {
	470594912, 650447616, 310934240, 695012864, 850358912,
61088076, 481306752, 786902080, 224042800, 805177664, 938284096,
145937936, 622867968, 160207584, 977329216, 716234240, 127727624,
415316352, 870137472, 18664444, 330872224, 93728752, 914779200,
736261248, 643647616, 755802688, 213052336, 410240448, 218974736,
109419280, 178026128, 689569664
};
ENDVERBATIM

FUNCTION fseed(seed) {
VERBATIM
    initstate((unsigned)_lseed,(char *)state2,32);
	_lfseed = _lseed;
ENDVERBATIM
}


FUNCTION n_rand() { : 0.0 <= n_rand < 1.0
VERBATIM
    _ln_rand = ((double)random()) / (((double)MAXLONG) + 1.);
ENDVERBATIM
}


FUNCTION fran(l, h) { : returns random number between low and high
VERBATIM
{
	int low, high;
    double num, imax, *getarg();
    
	low = (int)_ll;
	high = (int)_lh;
    imax = high-low+1; /* the total number of numbers being used */
    _lfran = (double)(low + (int) (imax*n_rand()));  
}
ENDVERBATIM
}


FUNCTION u_rand() { : uniform distribution (0.0 <= rand <= 1.0)
VERBATIM
    _lu_rand = (((double)random()) / ((double)MAXLONG));
ENDVERBATIM
}
    

FUNCTION norm() { : gaussian distribution around 0
VERBATIM
{
    static int iset = 0;
    static float gset;
    float fac, r , v1, v2;
    double sqrt();

    if (iset == 0) {
        do {
	    	v1 = 2.0 * n_rand() - 1.0;
		    v2 = 2.0 * n_rand() - 1.0;
		    r = v1 * v1 + v2 * v2;
	    } while (r >= 1.0);

        fac = (float)sqrt(-2.0 * log(r) / r);
        gset = v1 * fac;
        iset = 1;
        _lnorm = v2 * fac;

    } else {
        iset = 0;
        _lnorm = (double)gset;
    }
}
ENDVERBATIM
}


FUNCTION pois(mean) { : poisson distribution
VERBATIM
    _lpois = - _lmean * log(((double)random()+1.) / ((double)MAXLONG+1.));
ENDVERBATIM
}

FUNCTION poisint(mean) {
  poisint = poisrand(mean)
}
