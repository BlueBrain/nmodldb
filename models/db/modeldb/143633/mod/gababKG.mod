TITLE minimal model of GABAb receptors

COMMENT
-----------------------------------------------------------------------------

M. Birdno added kinetic scheme for effects due to Adenosine &/or Muscarinic ACh March 2009.
See Steriade, Jones, and McCormick (1997) Thalamus, Volume I (pp. 727-730)


	Kinetic model of GABA-B receptors
	=================================

  MODEL OF SECOND-ORDER G-PROTEIN TRANSDUCTION AND FAST K+ OPENING
  WITH COOPERATIVITY OF G-PROTEIN BINDING TO K+ CHANNEL

  PULSE OF TRANSMITTER

  SIMPLE KINETICS WITH NO DESENSITIZATION

	Features:

	  - peak at 100 ms; time course fit to Tom Otis' PSC
	  - SUMMATION (psc is much stronger with bursts)


	Approximations:

	  - single binding site on receptor	
	  - model of alpha G-protein activation (direct) of K+ channel
	  - G-protein dynamics is second-order; simplified as follows:
		- saturating receptor
		- no desensitization
		- Michaelis-Menten of receptor for G-protein production
		- "resting" G-protein is in excess
		- Quasi-stat of intermediate enzymatic forms
	  - binding on K+ channel is fast


	Kinetic Equations:

	  dR/dt = K1 * T * (1-R-D) - K2 * R

	  dS/dt = K5 * T * (1 - S) - K6 * S  

	  dG/dt = K3 * (R+S) - K4 * G

	  S : activated receptor due to neuromodulation by mACh &/or adenosine 
	  R : activated receptor
	  T : transmitter
	  G : activated G-protein
	  K1,K2,K3,K4,K5,K6 = kinetic rate cst

  n activated G-protein bind to a K+ channel:

	n G + C <-> O		(Alpha,Beta)

  If the binding is fast, the fraction of open channels is given by:

	O = G^n / ( G^n + KD )

  where KD = Beta / Alpha is the dissociation constant

-----------------------------------------------------------------------------

  Parameters estimated from patch clamp recordings of GABAB PSP's in
  rat hippocampal slices (Otis et al, J. Physiol. 463: 391-407, 1993).

-----------------------------------------------------------------------------

  PULSE MECHANISM

  Kinetic synapse with release mechanism as a pulse.  

  Warning: for this mechanism to be equivalent to the model with diffusion 
  of transmitter, small pulses must be used...

  For a detailed model of GABAB:

  Destexhe, A. and Sejnowski, T.J.  G-protein activation kinetics and
  spill-over of GABA may account for differences between inhibitory responses
  in the hippocampus and thalamus.  Proc. Natl. Acad. Sci. USA  92:
  9515-9519, 1995.

  For a review of models of synaptic currents:

  Destexhe, A., Mainen, Z.F. and Sejnowski, T.J.  Kinetic models of 
  synaptic transmission.  In: Methods in Neuronal Modeling (2nd edition; 
  edited by Koch, C. and Segev, I.), MIT press, Cambridge, 1996.

  This simplified model was introduced in:

  Destexhe, A., Bal, T., McCormick, D.A. and Sejnowski, T.J.
  Ionic mechanisms underlying synchronized oscillations and propagating
  waves in a model of ferret thalamic slices. Journal of Neurophysiology
  76: 2049-2070, 1996.  

  See also http://cns.iaf.cnrs-gif.fr

  Alain Destexhe, Salk Institute and Laval University, 1995
  27-11-2002: the pulse is implemented using a counter, which is more
       stable numerically (thanks to Yann LeFranc)

-----------------------------------------------------------------------------
ENDCOMMENT



INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS GABAbKG
	POINTER pre, vext, pmodyn
	USEION k READ ek WRITE ik
	RANGE C, R, S, G, g, gmax, lastrelease, TimeCount, ek, K1, K2, K5, K6, nsm, i
	GLOBAL Cmax, Cdur, Prethresh, Deadtime
	GLOBAL K3, K4, KD
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
	Prethresh = 0 			: voltage level nec for release
	Deadtime = 1	(ms)		: mimimum time between release events
	stimon = 0
	nsm = 1
:
:	From Kfit with long pulse (5ms 0.5mM)
:
	K1	= 0.52   (/ms mM) : 	: forward binding rate to receptor
	K2	= 0.0013 (/ms)	: 	: backward (unbinding) rate of receptor
	K3	= 0.098 (/ms)		: 	: rate of G-protein production
	K4	= 0.033 (/ms)		: 	: rate of G-protein decay
	K5	= 0.52   (/ms) : 	: forward binding rate to adenosine/mACh receptor
	K6	= 0.00013 (/ms)	: 	: backward (unbinding) rate of adenosine/mACh receptor
	KD	= 100			: dissociation constant of K+ channel
	n	= 4			: nb of binding sites of G-protein on K+
	gmax		(umho)		: maximum conductance
}

ASSIGNED {
	ek 		(mV)
	v		(mV)		: postsynaptic voltage
	ik 		(nA)		: current = g*(v - ek)
	g 		(umho)		: conductance
	C		(mM)		: transmitter concentration
	Gn
	pre 				: pointer to presynaptic variable
	lastrelease	(ms)		: time of last spike
	TimeCount	(ms)		: time counter
	vext				: vext turns into a dummy variable to determine whether a stim pulse is ON
	pmodyn
	i 		(nA)
}


STATE {
	R				: fraction of activated GABAb receptor
	S				: fraction of activated adenosine/mACh receptor evoked by STIMULATION
	G				: fraction of activated G-protein
}


INITIAL {
	C = 0
	lastrelease = -1000

	R = 0
	S = 0
	G = 0
	TimeCount=-1
}

BREAKPOINT {
	SOLVE bindkin METHOD euler
	Gn = G^n
	g = gmax * Gn / (Gn+KD)
	i = g*(v - ek) 
	ik = i
}


DERIVATIVE bindkin {

	release()		: evaluate the variable C

	R' = K1 * C * (1-R) - K2 * R
	S' = K5 * pmodyn * (1-S) - K6 * S
	G' = K3 * (R + 2 * nsm * S) - K4 * G
}


PROCEDURE release() {
	:will crash if user hasn't set pre with the connect statement 

	TimeCount=TimeCount-dt			: time since last release ended

						: ready for another release?
	if (TimeCount < -Deadtime) {
		if (pre > Prethresh) {		: spike occured?
			C = Cmax			: start new release
			lastrelease = t
			TimeCount=Cdur
		}
						
	} else if (TimeCount > 0) {		: still releasing?
	
		: do nothing
	
	} else if (C == Cmax) {			: in dead time after release
		C = 0.
	}
	
	if(vext==0) {      : vext turns into a dummy variable to determine whether a stim pulse is ON
		stimon=0
	}else{
		stimon=1
	}


}
