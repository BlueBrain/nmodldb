COMMENT
  ================================================================================
   Random number generator routines
 
   from C recipes Chapter 7. pp 204 ff.
   and Knuth, Seminumerical algs
   Adapted by Jose Ambros-Ingerson jose@kiubo.net

  ================================================================================
ENDCOMMENT

NEURON {
  SUFFIX nothing
}

VERBATIM
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h> /* contains MAXLONG */

/* Michael Hines fix for cygwin on mswin */
#if !defined(MAXLONG)
#include <limits.h>
#define MAXLONG LONG_MAX
#endif

static int  dev0_first = 1;
static long dev0_seed;

extern double du_dev0();
extern double du_dev( double min, double max );
extern int iu_dev( int min, int max );
extern double dexp_dev0();
extern double dexp_dev( double lambda );
extern double dgauss_dev0();
extern double dgauss_dev( double mean, double sdev );
extern int igeom_dev( double pg );
double ran2();	/* ran2 not to be used elsewhere */

ENDVERBATIM

: ================================================================================
: set dev0_seed; we assume seed is positive. If not ran2 will check and correct
: make accesible from hoc 
FUNCTION set_dev0_seed( seed ){
VERBATIM
  dev0_first = 0;
  dev0_seed = -1 * (long) _lseed;
  ran2( &dev0_seed );
ENDVERBATIM
  set_dev0_seed = 1
}

VERBATIM
/* ================================================================================
   double uniform deviate in (0,1)
 ================================================================================ */
double du_dev0(){
  double r2;
  if( dev0_first ){
    dev0_first = 0;
    dev0_seed = 1234;
  }
  r2 = ran2( &dev0_seed );
  return r2;
}
ENDVERBATIM

VERBATIM
/* ================================================================================
   return a uniform deviate in (min,max)
   ================================================================================ */
double du_dev( double min, double max )
{
  return( min + (max-min)*du_dev0() );
}

/* ================================================================================
   return an integer uniform deviate in [min,max]
   ================================================================================ */
int iu_dev( int min, int max )
{
  return( (int)((double)min + (double)(max-min+1)*du_dev0()) );
}

/* ================================================================================
   Returns an exponentially distributed, positive random deviate of unit mean,
   using du_dev1 as the source of uniform deviates
   ================================================================================ */
double dexp_dev0(){
  double dum; 
  do dum=du_dev0(); while ( dum==0.0);
  return( -log( dum ));
}

/* ================================================================================
   p(x) = lambda e^{-lambda x}
   ================================================================================ */
double dexp_dev( double lambda )
{
  return( dexp_dev0()/lambda );
}

/*======================================================================
  generate pseudorandom numbers according to normal distribution with
  mean 0 and std dev 1.
  We use algorithm C, page 117 of Knuth's Seminumerical Algorithms
  ======================================================================*/
double dgauss_dev0()
{
  static int    NotHave=1;
  static double X2;
  double        V1, V2, S, LS;

  if( NotHave ) {
    NotHave = 0;
    do {
      V1 = 2.0 * du_dev0() - 1.0;
      V2 = 2.0 * du_dev0() - 1.0;
    } while( (S = V1*V1 + V2*V2) > 1.0 );

    LS = sqrt( (-2.0 * log( S )) / S );
    X2 = V2 * LS;
    return( V1 * LS );
  }
  else {
    NotHave = 1;
    return( X2 );
  }
}

/*======================================================================
  generate pseudorandom numbers according to normal distribution with
  mean mean and std dev sdev.
  We use algorithm C, page 117 of Knuth's Seminumerical Algorithms
  ======================================================================*/
double dgauss_dev( double mean, double sdev )
{
  return( dgauss_dev0()*sdev + mean );
}

/*======================================================================
  generate pseudorandom numbers according to a geometic distribution
  with prob 0<pg<=1
  ======================================================================*/
int igeom_dev( double pg ){
  int ig;
  ig = 1;
  while( pg < du_dev0() && pg<=1.0 ) ig++;
  return( ig );
}

/* ================================================================================
   From Numerical Recipes 
  Long period (2> 2 x 10^18) random number generator of L'Ecuyer with Bays-Durham shuffle
  and added safeguards. Returns a uniform deviate in (0.0,1.0). 
    Call with idum a negative integer to initialize; thereafter, do not alter idum
  between succesive deviates in a sequence. 
    RNMX should approximate the largest floating point that is less than 1
 ================================================================================ */
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  double temp;
  
  if (*idum <= 0) {			/* Initialize */
    if (-(*idum) < 1) *idum=1;		/* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) {		/* Load shuffle table (after 8 warm-ups) */
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;			/* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
ENDVERBATIM
