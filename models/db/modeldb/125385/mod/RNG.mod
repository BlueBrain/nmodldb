TITLE RNG.mod

COMMENT
----------------------------------------------------------------
Random Number Generator based on Numerical Recipes implementation
of ran1(), gasdev(), and bnldev().

Functions UniDev_RNG, BinDev_RNG, and GasDev_RNG provide deviates
drawn from this distribution.

The seed to be used at the beginning of the simulation is set by
giving the parameter init_seed the desired value. If set to -1, then
the generator continues to run across the initialization, essentially
starting with a random seed.

SetSeedNow_RNG sets the initial seed to the given value.

******** WARNING ********
There can only be one of these point processes in any given simulation.
This is due to the use of static arrays of integers, which can not
easily be allocated locally within each segment.

This file must have the same vectorization state as all other 
models it is linked with so that the functions and procedures
are callable from other methods.

Additional warning. This code uses knowledge of specific name-munging 
performed by neuron when passing arguments to underlying C functions.
This is highly non-portable. Neuron version 3 labels arguments and
local variables, named XX in MOD code as _lXX in the generated C code.
Each occurence is commented with **NP1** for ease of location.
*************************

Note several modifications from the original Numerical Recipes source:
1. idum is a static variable of uniran1() and is not sent as an argument.
2. gasdev() and bnldev() does not pass idum to ran1().
3. Since bnldev have underflow errors when computing exp(oldg-gammln(em+1.0)
      -gammln(en-em+1.0)+em*plog+(en-em)*pclog)  
    this set errno, and causes neuron to print an error message.
    Thus errno is checked before and after this computation and if 
    changed it being set back to its old value.

Jan 1999 Miki London & Peter N. Steinmetz
May 2002 Kamran Diba changed init_seed to RANGE variable.
August 2002 all the errno business was removed to work with Neuron 5.2 for windows
----------------------------------------------------------------
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
    init_seed = -1           : seed to use when simulation initialized
                        :  unless set to -1, which lets the RNG
                        :  free run across simulations
    DONT_VECTORIZE          : required declaration
}

NEURON {
    POINT_PROCESS RNG
    RANGE init_seed
    GLOBAL DONT_VECTORIZE   : inhibits vectorization
}

INITIAL {
VERBATIM
    /* use init_seed, unless is set to -1 allowing free running */
    if (init_seed != -1)
        set_seed(init_seed);
    /* otherwise leave it alone */
ENDVERBATIM
}

VERBATIM
/* ----------------------------------------------------------------
   Block of definitions and functions to implement Numerical Recipes
   generators.
   ---------------------------------------------------------------- */
#include <math.h>
#include <limits.h>

#define        PI 3.141592654
#define        r_ia     16807
#define        r_im     2147483647
#define        r_am     (1.0/r_im)
#define        r_iq     127773
#define        r_ir     2836
#define        r_ntab   32 
#define        r_ndiv   (1+(r_im-1)/r_ntab)
#define        r_eps    1.2e-7
#define        r_rnmx   (1.0-r_eps)

/* ---------------------------------------------------------------- */
/* gammln - compute natural log of gamma function of xx */
static double
gammln(double xx)
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

/* ---------------------------------------------------------------- */
/* uniran1 static globals */
static long idum_RNG = -1;
static long iy_RNG = 0; 
static long iv_RNG[r_ntab];

/* ---------------------------------------------------------------- */
/* set_seed - setup seed and temporaries for uniran rng */
static void
set_seed(long seed) {
    int j, k;
    
    idum_RNG = fabs(seed);
    fprintf(stderr,"RNG: seed is set to %ld\n",idum_RNG); 

    for ( j = r_ntab + 7; j>=0;j--){
        k = (idum_RNG/r_iq);
        idum_RNG = r_ia*(idum_RNG-k*r_iq)-r_ir*k;
        if (idum_RNG<0) idum_RNG += r_im;
        if (j<r_ntab) iv_RNG[j] = idum_RNG;
    }
    iy_RNG = iv_RNG[0];
}

/* ---------------------------------------------------------------- */
/* uniran1 - draw a random deviate uniformly distributed on 0 - 1

   This is a modification of ran1() to be stand alone with its own 
   static global variables. 
   Call set_seed to setup the seed now.
 */

static double
uniran1()
{     
    int  j;
     long k;
     double temp; 
    k = (idum_RNG/r_iq);
     idum_RNG = r_ia*(idum_RNG-k*r_iq)-r_ir*k;
     if (idum_RNG<0) idum_RNG += r_im;
     j = iy_RNG/r_ndiv;
     iy_RNG = iv_RNG[j];
     iv_RNG[j] = idum_RNG;
     if ((temp = r_am*iy_RNG) >r_rnmx) return r_rnmx;
     else return temp;
}

/* ---------------------------------------------------------------- */
/* gasdev - draw gaussian random deviate */
static double
gasdev() {
     double fac,r,v1,v2;
    static int iset_RNG = 0;
    static double gset_RNG;

     if  (iset_RNG == 0) {
        do {
            v1=2.0*uniran1()-1.0;
               v2=2.0*uniran1()-1.0;
               r=v1*v1+v2*v2;
        } while (r >= 1.0 || r == 0.0);
        fac=sqrt(-2.0*log(r)/r);
        gset_RNG=v1*fac;
        iset_RNG=1;
        return v2*fac;
    } else {
        iset_RNG=0;
          return gset_RNG;
    }
}

/* ---------------------------------------------------------------- */
/* bnldev - draw binomial random deviate */
double
bnldev(double ppr, int nnr) {
    int j;
    static int nold=(-1);
    double am,em,g,angle,p,bnl,sq,bt,y;
    static double pold=(-1.0),pc,plog,pclog,en,oldg;
    
    /* prepare to always ignore errors within this routine */
     
    
    p=(ppr <= 0.5 ? ppr : 1.0-ppr);
    am=nnr*p;
    if (nnr < 25) {
        bnl=0.0;
        for (j=1;j<=nnr;j++)
            if (uniran1() < p) bnl += 1.0;
    }
    else if (am < 1.0) {
        g=exp(-am);
        bt=1.0;
        for (j=0;j<=nnr;j++) {
            bt *= uniran1();
            if (bt < g) break;
        }
        bnl=(j <= nnr ? j : nnr);
    }
    else {
        if (nnr != nold) {
            en=nnr;
            oldg=gammln(en+1.0);
            nold=nnr;
        }
        if (p != pold) {
            pc=1.0-p;
             plog=log(p);
            pclog=log(pc);
            pold=p;
        }
        sq=sqrt(2.0*am*pc);
        do {
            do {
                angle=PI*uniran1();
                    angle=PI*uniran1();
                y=tan(angle);
                em=sq*y+am;
            } while (em < 0.0 || em >= (en+1.0));
            em=floor(em);
                bt=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0) - 
                gammln(en-em+1.0)+em*plog+(en-em)*pclog);
        } while (uniran1() > bt);
        bnl=em;
    }
    if (p != ppr) bnl=nnr-bnl;
    
    /* recover error if changed during this routine, thus ignoring
        any errors during this routine */
   
    
    return bnl;
}

#undef PI 
#undef r_ia 
#undef r_im 
#undef r_am 
#undef r_iq  
#undef r_ir 
#undef r_ntab 
#undef r_ndiv 
#undef r_eps 
#undef r_rnmx 
ENDVERBATIM

: ----------------------------------------------------------------
: SetSeedNow - set seed for generator immediately
PROCEDURE SetSeedNow(seed) {
VERBATIM
    set_seed((long)_lseed); /* **NP1** non-portable argument use */
ENDVERBATIM
}

: ----------------------------------------------------------------
: UniDev - draw a uniform deviate from the generator
FUNCTION UniDev () {
VERBATIM
    _lUniDev = uniran1();   /*  **NP1** non-portable argument use */
ENDVERBATIM
}
    
: ----------------------------------------------------------------
: BnlDev - draw a uniform deviate from the generator
FUNCTION BnlDev (ppr, nnr) {
VERBATIM
    _lBnlDev = bnldev(_lppr, (int)_lnnr);  /*  **NP1**
                                         non-portable argument use */
ENDVERBATIM
}

: ----------------------------------------------------------------
: test generator
PROCEDURE test_ran1() {
    LOCAL j,x
      VERBATIM 
      FILE *out; 
      out = fopen("test_ran1.out","w");
      ENDVERBATIM 

    j = 0
    while (j<10000) {
      VERBATIM 
       fprintf(out,"%lg\n",uniran1());
      ENDVERBATIM 
      j = j + 1
    }
      VERBATIM 
      fclose(out);
      ENDVERBATIM 
   
}       

: ----------------------------------------------------------------
: test_gasdev - test gaussian deviates generator
PROCEDURE test_gasdev() {
    LOCAL j,x
      VERBATIM 
      FILE *out; 
      out = fopen("test_gasdev.out","w");
      ENDVERBATIM 

    j = 0
    while (j<10000) {
      VERBATIM 
       fprintf(out,"%lg\n",gasdev());
      ENDVERBATIM 
      j = j + 1
    }
      VERBATIM 
      fclose(out);
      ENDVERBATIM 
   
}       
