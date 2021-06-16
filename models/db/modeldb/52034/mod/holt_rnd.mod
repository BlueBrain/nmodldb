TITLE Random number generator
COMMENT
$Id: holt_rnd.mod,v 1.1.1.2 1993/08/07 15:36:29 holt Exp $

Revision 1.1  1993/08/04  14:36:28  holt
Initial revision

These functions generate random numbers with various distributions.
They all call the function holt_random(), which is a random number generator
which decides which of the standard random number generators to call.
The seed value is the time if not explicitly specified.
 
ENDCOMMENT

NEURON {
        SUFFIX random
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

VERBATIM
double holt_random();
void holt_seed();
double holt_normrand();
double holt_exprand();
int holt_poisrand();
ENDVERBATIM

:
: Returns uniformly distributed random numbers on [0,1).
:
FUNCTION Uniform() {
  Uniform = holt_random()	: Return a random number.
}

:
: Returns Poisson-distributed integers from a Poisson distribution with
: the specified mean.
:
FUNCTION Poisson(mean) {
VERBATIM
  _lPoisson = holt_poisrand(_lmean);
ENDVERBATIM
}

:
: Returns a random number from an exponential distribution.
:
FUNCTION Exp() {
  Exp = holt_exprand()		: Return a random number.
}

:
: Returns a random number from a normal distribution.
:
FUNCTION Normal(mean, std_dev) {
VERBATIM
  _lNormal = holt_normrand(_lmean, _lstd_dev); /* Must do in VERBATIM block or
					        * else another argument gets
					        * added in front. */
ENDVERBATIM
}

:
: Seed the random number generator, and specify which algorithm to use:
: Arguments:
: 1) The seed value (an integer).
: 2) The algorithm.  0 = drand48(), 1 = random().
:
PROCEDURE Seed(seedval, algorithm) {
VERBATIM
  holt_seed((int)_lseedval, (int)_lalgorithm);
ENDVERBATIM
}

VERBATIM
/* A collection of routines for generating random numbers.
 * I abstracted some of these routines from the standard scopmath library
 * from Duke university so they could use a different random number generator.
 */

int seed = 0;			/* Seed for the random number generator. */
#define DRAND48 0		/* Which random number generator to use. */
#ifndef hpux			/* HP/UX doesn't come with random(). */
#define RANDOM 1
#define N_GENERATORS 2		/* How many different generators are implemented. */
#else
#define N_GENERATORS 1
#endif
#define DEFAULT_GENERATOR DRAND48 /* Which is used by default. */
int random_generator;		/* Code for which generator to use. */

#ifdef hpux
extern double drand48();
#endif
#ifndef __alpha
extern long random();
extern double drand48();
#endif


/*
 * Function to seed the random number generator:
 * Arguments:
 * 1) The seed.  If this is 0, the time is used for the seed.
 * 2) Which random number generator to use.
 *    0 = drand48()
 *    1 = random()
 */
void
holt_seed(seedval, gen_code)
int seedval, gen_code;
{
  if (seedval == 0)		/* No seed explicitly given? */
	seedval=3491;
//    seedval = time(0);		/* Seed based on the time. uncomment on unix/linux systems */

  seed = seedval;

  if (gen_code >= N_GENERATORS)	/* Illegal parameter? */
    gen_code = DEFAULT_GENERATOR; /* Don't give an error message. */
  random_generator = gen_code;	/* Remember which generator to use. */

  switch (random_generator)	/* Seed is different depending on which */
  {				/* algorithm we use. */
  case DRAND48:    srand48((long)seedval);    break;
// #ifndef hpux			// can uncomment these three lines on linux/unix systems if desired
//  case RANDOM:     srandom(seedval);          break; 
// #endif
  }
}

/*
 * Function to return a uniformly distributed random number on [0,1).
 * Seeds itself automatically with the time if no seed has been provided.
 */
double
holt_random()
{
  if (seed == 0)		/* No seed provided yet? */
    holt_seed(0, DEFAULT_GENERATOR); /* Seed the random number generator. */

  switch (random_generator)	/* Use the proper generator: */
  {
  case DRAND48:   return drand48();
// #ifndef hpux
//  case RANDOM:    return (double)random() / (double)0x7fffffff;
// #endif
  }
}

/*--------------------------------------------------------------*/
/*								*/
/*  NORMRAND							*/
/*								*/
/*    Selects a random number from the normal distribution with	*/
/*    specified mean and standard deviation			*/
/*								*/
/*  Returns: Double precision number from the distribution	*/
/*								*/
/*  Calling sequence: normrand(mean, std_dev)			*/
/*								*/
/*  Arguments							*/
/*    Input:	mean, double	mean of the distribution	*/
/*								*/
/*		std_dev, double	standard deviation of the	*/
/*				distribution			*/
/*								*/
/*    Output:	arguments are unchanged				*/
/*								*/
/*  Functions called: random					*/
/*								*/
/*--------------------------------------------------------------*/

double 
holt_normrand(mean, std_dev)
double mean, std_dev;
{
    double s, v1, v2;
    extern double sqrt(), log();

    s = 1.0;
    while (s >= 1.0)
    {
	v1 = 2.0 * holt_random() - 1.0;
	v2 = 2.0 * holt_random() - 1.0;
	s = (v1 * v1) + (v2 * v2);
    }
    v2 = v1 * sqrt(-2.0 * log(s) / s);
    return (v2 * std_dev + mean);
}

/*
 * This function returns the number of events in a poisson process
 * with mean xm during one iteration.  Taken from numerical recipes
 * in C.
 */
#define PI 3.141592654

int
holt_poisrand(xm)
double xm;
{
  float gammln();
  double drand48();
  static float sq,alxm,g,oldm=(-1.0);
  double em,tt,y;
  double exp();

  if (xm < 12.0) {
    if (xm != oldm) {
      oldm=xm;
      g=exp(-xm);
    }
    em = -1;
    tt=1.0;
    do {
      ++em;
      /*                      tt *= ran1(idum);  */
      tt *= holt_random();
    } while (tt > g);
  } else {
    if (xm != oldm) {
      oldm=xm;
      sq=sqrt(2.0*xm);
      alxm=log(xm);
      g=xm*alxm-gammln(xm+1.0);
    }
    do {
      do {
	/*              y=tan(PI*ran1(idum));  */
	y=tan(PI*holt_random());
	em=sq*y+xm;
      } while (em < 0.0);
      em=floor(em);
      tt=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
    } while (holt_random() > tt);
  }
  return (int)em;
}
#undef PI
/* (C) Copr. 1986-92 Numerical Recipes Software #.3. */

float
gammln(xx)
float xx;
{
        double x,tmp,ser;
        static double cof[6]={76.18009173,-86.50532033,24.01409822,
                -1.231739516,0.120858003e-2,-0.536382e-5};
        int j;

        x=xx-1.0;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.0;
        for (j=0;j<=5;j++) {
                x += 1.0;
                ser += cof[j]/x;
        }
        return -tmp+log(2.50662827465*ser);
}

/****************************************************************/
/*								*/
/*  Abstract: exprand()						*/
/*								*/
/*    Random sample from the exponential distribution exp(-x)	*/
/*								*/
/*  Calling sequence: exprand()					*/
/*								*/
/*  Arguments:							*/
/*	Input:	none						*/
/*								*/
/*	Output:	none						*/
/*								*/
/*  Functions called: holt_random()				*/
/*								*/
/*  Returns: Double precision number from the distribution	*/
/*								*/
/*  Files accessed: none					*/
/*								*/
/*								*/
/****************************************************************/

double 
holt_exprand()
{
    extern double log();

    return (-log(holt_random()));
}
ENDVERBATIM


