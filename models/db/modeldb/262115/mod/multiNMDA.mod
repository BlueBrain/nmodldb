TITLE multiple NMDA receptors

COMMENT
-----------------------------------------------------------------------------
Mechanism for handling multiple presynaptic entries to the same compartment;
up to 1000 synapses can be handled using a different pointer that must be set
for each presynaptic variable using addlink.  Optimization algorithm from
Lytton, W.W., Neural Computation, in press, 1996.  This mechanism allows
considerable acceleration of the simulation time if many receptors of the same
type must be simulated in the same compartment.

This file was configured for NMDA receptors.  The mechanism was a first-order
kinetic model with pulse of transmitter (see Destexhe, A., Mainen, Z. and
Sejnowski, T.J.  Neural Computation, 6: 14-18, 1994).

Parameters were obtained from fitting the model to whole-cell recorded NMDA
postsynaptic currents (Hessler et al., Nature 366: 569-572, 1993).  The fit
was performed using a simplex algorithm using short pulses of transmitter (0.5
mM during 0.3 ms).

-----------------------------------------------------------------------------
EXAMPLE OF HOW TO USE:

create POST,PRE[10]		// create compartments
objectvar c			// create an object
c = new multiNMDA()		// create multiple NMDA kinetic synapses
POST c.loc(0.5)			// localize synapse on postsyn compartment
c.gmax = 0.001			// assign max conductance of each syn (mu S)
c.allocate(10)			// allocate space for 10 presyn variables
for i=0,9 { 			// link presynaptic variables
   c.addlink(&PRE[i].v)
}  
-----------------------------------------------------------------------------
WARNINGS:

  - only ok for synaptic mechanisms where all weights are equal
    (see Lytton paper for implementation of different weights)


  Alain Destexhe, Laval University, 1995

-----------------------------------------------------------------------------
ENDCOMMENT

: defines maximal number of possible links to presynaptic variables
: this number should correpond to the number of pointers pre00, pre01, ...
: defined in the NEURON block

DEFINE MAXSYNNMDA 250
VERBATIM
static int MAXSYNNMDA = 250;
ENDVERBATIM

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS multiNMDA
	NONSPECIFIC_CURRENT i
	RANGE Ron, Roff, ri, nsyn, non, g, gmax, B
	GLOBAL Cmax, Cdur, Alpha, Beta, Erev
	GLOBAL Prethresh, Deadtime, Rinf, Rtau, mg
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {
	dt		(ms)

	Cmax	= 0.5	(mM)		: max transmitter concentration
	Cdur	= 0.3	(ms)		: transmitter duration (rising phase)
	Alpha	= 0.11	(/ms mM)	: forward (binding) rate
	Beta	= 0.0064 (/ms)		: backward (unbinding) rate
	Erev	= 0	(mV)		: reversal potential
	gmax		(umho)		: maximum conductance of each synapse

	Prethresh = 0 			: voltage level nec for release
	Deadtime = 1	(ms)		: mimimum time between release events
	mg	= 1    (mM)		: external magnesium concentration
}


ASSIGNED {
	on[MAXSYNNMDA]			: state of each synapse
	TL[MAXSYNNMDA]	(ms)		: time since last event for each synapse
	ri[MAXSYNNMDA]			: state variable of each synapse
	lastrelease[MAXSYNNMDA] (ms)	: last release for each synapse
	Ron				: sum of all "on" synapses
	Roff				: sum of all "off" synapses
	nsyn				: number of synapses
	non				: number of synapses on

	Rinf				: steady state channels open
	Rtau		(ms)		: time constant of channel binding

	v		(mV)		: postsynaptic voltage
	i 		(nA)		: total current = g*(v - Erev)
	g 		(umho)		: total conductance

	trel		(ms)		: temp var
	ptr_array_nmda			: pointer array
	B				: magnesium block
}

INITIAL { LOCAL j
	FROM j=0 TO nsyn-1 {
		on[j] = 0
		TL[j] = -9e9
		lastrelease[j] = -9e9
		ri[j] = 0
	}
	Ron = 0
	Roff = 0
	non = 0

	Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)
	Rtau = 1 / ((Alpha * Cmax) + Beta)
}

BREAKPOINT {
   if(gmax > 0) {
	SOLVE release
	B = mgblock(v)		: B is the block by magnesium at this voltage
	g = gmax * (Ron+Roff) * B
	i = g*(v - Erev)
   } else {
	i = 0
   }
}

PROCEDURE release() { LOCAL q,j

  FROM j=0 TO nsyn-1 {	: update Ron, Roff, non for each synapse

    trel = ((t - lastrelease[j]) - Cdur)	: time since last release ended

    if (trel > Deadtime) {			: ready for another release?
				
	if (presynaptic(j) > Prethresh) {	: spike occured?
	  on[j] = 1			: start new release
	  non = non + 1
	  lastrelease[j] = t		: memorize release time
	  ri[j] = ri[j] * exptable( - Beta * (t-TL[j]))
					: evaluate state variable
	  TL[j] = t			: memorize last event
	  Ron = Ron + ri[j]		: increase Ron
	  Roff = Roff - ri[j]		: decrease Roff
	  if(Roff < 1e-9) { Roff = 0 }	: prevent roundoff errors
	}
						
    } else if (trel < 0) {			: still releasing?

		: do nothing
	
    } else if (on[j] > 0) {			: end of release ?
	on[j] = 0			: stop release
	non = non - 1
	ri[j] = Rinf + (ri[j]-Rinf) * exptable(- (t-TL[j]) / Rtau)
					: evaluate state variable
	TL[j] = t			: memorize last event
	Ron = Ron - ri[j]		: decrease Ron
	Roff = Roff + ri[j]		: increase Roff
	if(Ron < 1e-9) { Ron = 0 }		: prevent roundoff errors
    }

  }


  if(Roff > 0) {			: update Roff
     Roff = Roff * exptable(- Beta * dt)
     if(Roff < 1e-9) { Roff = 0 }	: prevent roundoff errors
  }

  if(non > 0) {				: update Ron
    q = non * Rinf
    Ron = q + (Ron - q) * exptable(- dt / Rtau) 
  }

}


FUNCTION exptable(x) { 
	TABLE  FROM -25 TO 25 WITH 10000

	if ((x > -25) && (x < 25)) {
		exptable = exp(x)
	} else {
		exptable = 0.
	}
}



:FUNCTION exptable(x) { 
:	TABLE  FROM -10 TO 10 WITH 10000
:
:	if ((x > -10) && (x < 10)) {
:		exptable = exp(x)
:	} else {
:		exptable = 0.
:	}
:}



:FUNCTION exptable(x) {
:	if(x > -50) {
:		exptable = exp(x)
:	} else {
:		exptable = 0.
:	}
:}



FUNCTION mgblock(v(mV)) {
	TABLE 
	DEPEND mg
	FROM -140 TO 80 WITH 1000

	mgblock = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}




:-------------------------------------------------------------------
:  Procedures for pointer arrays in nmodl 
:  create a pointer array and link its pointers to variables passed
:  from hoc (adapted from Mike Hines)
:-------------------------------------------------------------------


VERBATIM
#define ppnmda ((double***)(&(ptr_array_nmda)))
extern double* hoc_pgetarg();
ENDVERBATIM


:
: Procedure to allocate space for n pointers
:
PROCEDURE allocate(n) {
  VERBATIM
	if (*ppnmda) {
	   free(*ppnmda);
	}
	*ppnmda = ((double**) hoc_Ecalloc((int)_ln, sizeof(double *))), hoc_malchk();
  ENDVERBATIM
}

:
: procedure to get the value of a presynaptic variable
: index is the number of the presynaptic var
:
FUNCTION presynaptic(index) {
  VERBATIM
	if(_lindex >= nsyn) {
	   printf("Warning: attempt to use pointer outside range\n");
	   printf(" trying to use pointer number %d\n",(int)_lindex);
	   printf(" but number of defined pointers was nsyn=%d.\n",(int) nsyn);
	}
	_lpresynaptic = *((*ppnmda)[(int)_lindex]);
  ENDVERBATIM
}


:
: procedure to add a new presynaptic variable
: the address of the variable is passed as argument (from hoc)
: a new pointer is then linked to that variable
:
PROCEDURE addlink() {
  VERBATIM
	if(++nsyn > MAXSYNNMDA) {
	  printf("Exceeding maximum of allowed links MAXSYNNMDA=%d\n",MAXSYNNMDA);
	  printf("  edit the nmodl code to increase the maximum allowed.\n");
	  exit(-1);
	}
	(*ppnmda)[(int)(nsyn-1)] = hoc_pgetarg(1);
  ENDVERBATIM
}
