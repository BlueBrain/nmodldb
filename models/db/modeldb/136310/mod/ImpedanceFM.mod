TITLE Frequency-dependent impedance

COMMENT
----------------------------------------------------------------------

This mod files calculates the impedance of the extracellular medium
for the whole frequency range.  

There is no parameter to pass to the MOD file; but the procedure 
"calc_impedances" is callable from HOC.  The parameters of this
procedure are:
- fmax = max frequency at which the impedance must be calculated
- df = frequency step to calculate impedances
- rext = distance between the source and the electrode
- rmax = max distance to integrate
- dr = integration step for distance
- R = diameter of the current source
- sigma1 = extracellular conductivity close to the source (high)
	   (normalized value - normally equal to 1)
- sigma2 = "mean" extracellular conductivity far away (low)
	   (normalized value, fraction of the conductivity at the source)
- lambda = space constant of the exponential decay of conductivity
- epsilon = permittivity (normalized)
- sigmaR = absolute value of conductivity in S/um

(all distances are in um)

The resulting impedances Z are written in the 2 vector vec1 and
vec2, with vec1=real part of Z, vec2=imaginary part of Z.  These 
impedances constitute the "filter" of the extracellular space in
this model.

See details in:

  Bedard C, Kroger H & Destexhe A.  Modeling extracellular field 
  potentials and the frequency-filtering properties of extracellular 
  space.  Biophysical Journal 86: 1829-1842, 2004.

A. Destexhe, CNRS
destexhe@iaf.cnrs-gif.fr
see http://cns.iaf.cnrs-gif.fr

----------------------------------------------------------------------
ENDCOMMENT


NEURON	{ 
	POINT_PROCESS ImpedanceFM
        POINTER vec1
        POINTER vec2
}


ASSIGNED {
	vec1
	vec2
}


VERBATIM
#include <malloc.h>

// calculate impedance for the whole range of frequencies

void calc_impedances(Zr,Zi,fmax,df,rext,rmax,dr,R,sigma1,sigma2,lambda,epsilon,sigmaR)
	double Zr[],Zi[],fmax,df,rext,rmax,dr,R,sigma1,sigma2,lambda,epsilon,sigmaR;
	// vector Z is impedance, Z[0]=real part, Z[1]=imaginary part
   {
        double epsR,sigR,w,w2,sig,eps,den,ReZ,ImZ,r,f;
	float *sigmatab;
	double PI = 3.1415927;
	int j,k,siz;


	// printf("%s\n%g %g %g %g %g %g %g %g %g %g %g\n",\
	// "fmax,df,rext,rmax,dr,R,sigma1,sigma2,lambda,epsilon,sigmaR = ",\
	// fmax,df,rext,rmax,dr,R,sigma1,sigma2,lambda,epsilon,sigmaR);

	// tabulate all values of sigma in a table "sigmatab",
	// which avoids calculating them several times
	siz = fmax/df + 1;
	sigmatab = (float *) malloc(sizeof(float) * siz);
	k = 0;
	for(r=rext; r<=rmax; r=r+dr) {
	  sigmatab[k] = sigma2 + (sigma1-sigma2) * exp(-(r-R)/lambda);
	  k++;
	}

	// calculate the impedances for each frequency
	sigR = sigma1;
	epsR = epsilon;
	j=0;
	for(f=0; f<=fmax; f=f+df) {	// loop on frequencies
	  w = 2*PI*f;
	  w2 = w*w;
	  ReZ=0;
	  ImZ=0;
	  k=0;
	  for(r=rext; r<=rmax; r=r+dr) {	// compute integral
	    /* sig = sigmaF(r,R,sigma1,sigma2,lambda); */
	    sig = sigmatab[k];	// tabulated sigma
	    eps = epsilon;
	    den = r*r * (sig*sig + w2 * eps*eps);
	    ReZ = ReZ + (sig*sigR + w2*eps*epsR) / den;
	    ImZ = ImZ + (sig*epsR - sigR*eps) / den;
	    k++;
	  }
	  Zr[j] = dr/(4*PI*sigmaR) * ReZ;	// impedance (UNITS: Ohm)
	  Zi[j] = w * dr/(4*PI*sigmaR) * ImZ;
          j++;
	  /* printf("last sigma = %g\n",sig); */
	}
	free(sigmatab);
   }


ENDVERBATIM





:
: procedure callable from hoc, and which calls the mutual information
: algorithm.  The address of the pointed vecors are passed to the
: C procedure above.
:
PROCEDURE impedance(fmax,df,rext,rmax,dr,R,sigma1,sigma2,lambda,epsilon,sigmaR) {

VERBATIM

  /* printf("argument passed in procedure : f = %g\n",_lf); */

  calc_impedances(&vec1,&vec2,_lfmax,_ldf,_lrext,_lrmax,_ldr,_lR,_lsigma1,\
	_lsigma2,_llambda,_lepsilon,_lsigmaR);

ENDVERBATIM
}

