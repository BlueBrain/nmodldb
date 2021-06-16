: $Id: nmda.mod,v 1.4 1994/05/24 15:56:52 billl Exp $
COMMENT
-----------------------------------------------------------------------------
Simple synaptic mechanism derived for first order kinetics of
binding of transmitter to postsynaptic receptors.

A. Destexhe & Z. Mainen, The Salk Institute, March 12, 1993.

Last modif.  March 4th, 1994.


Reference:

Destexhe, A., Mainen, Z. and Sejnowski, T.J.  An efficient method for 
computing synaptic conductances based on a kinetic model of receptor binding.
Neural Computation, 6: 14-18, 1994.

-----------------------------------------------------------------------------
  Simple NMDA based on GABA SYNAPSE (GABA-A receptors) in gabalow.mod
-----------------------------------------------------------------------------
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  POINT_PROCESS NMDA
  POINTER pre
  RANGE C, R, R0, R1, g, gmax, lastrelease, B, spk
  NONSPECIFIC_CURRENT i
  GLOBAL Cmax, Cdur, Alpha, Beta, Erev, Prethresh, Deadtime, Rinf, Rtau
}

INCLUDE "queue.inc"  : queue routines

UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
  (umho) = (micromho)
  (mM) = (milli/liter)
}

PARAMETER {
  mg    = 1.    (mM)            : external magnesium concentration
  Cmax	= 1	(mM)		: max transmitter concentration
  Cdur	= 5.0	(ms)		: transmitter duration (rising phase)
  Alpha	= 0.072	(/ms mM)	: forward (binding) rate
  Beta	= 0.0066 (/ms)		: backward (unbinding) rate 1/150
  Erev	= 0	(mV)		: reversal potential
  Prethresh = 0 		: voltage level nec for release
  Deadtime = 1	(ms)		: mimimum time between release events
  gmax	= 0.001 (umho)		: maximum conductance
  vmin = -120     (mV)
  vmax = 100      (mV)
}


ASSIGNED {
  v		(mV)		: postsynaptic voltage
  i 		(nA)		: current = g*(v - Erev)
  g 		(umho)		: conductance
  C		(mM)		: transmitter concentration
  B                             : fraction free of Mg2+ block
  R				: fraction of open channels
  R0				: open channels at start of release
  R1				: open channels at end of release
  Rinf				: steady state channels open
  Rtau		(ms)		: time constant of channel binding
  pre 				: pointer to presynaptic variable
  spk                           : flag for spk occuring
  lastrelease	(ms)		: time of last spike
}

INITIAL {
  initq()                       : ****

  R = 0
  C = 0
  R0 = 0
  R1 = 0
  Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)
  Rtau = 1 / ((Alpha * Cmax) + Beta)
  lastrelease = -9e9
  rates(v)
}

BREAKPOINT {
  rates(v)
  SOLVE release
  g = gmax * R * B
  i = g*(v - Erev)
}

PROCEDURE release() { LOCAL q
  :will crash if user hasn't set pre with the connect statement 

  if (! spk && pre > Prethresh) { : new spike occured?
    spk = 1
    pushq(t+delay) } 

  if (spk && pre < Prethresh) { : spike over?
    spk = 0 }

  q = ((t - lastrelease) - Cdur) : time since last release ended

  : ready for another release?
  if (q > Deadtime) {

    if (t >= queu[head]) {      : **** a current spike time
      popq()                    : ****
      C = Cmax			: start new release
      R0 = R
      lastrelease = t
    }
    
  } else if (q < 0) {		: still releasing?

    if (t > queu[head]) { popq() } : **** throw away value from missed spikes

  } else if (C == Cmax) {	: in dead time after release
    R1 = R
    C = 0.
  }
  
  if (C > 0) {			: transmitter being released?
    R = Rinf + (R0 - Rinf) * exptable (- (t - lastrelease) / Rtau)
  } else {			: no release occuring
    R = R1 * exptable (- Beta * (t - (lastrelease + Cdur)))
  }
  
  VERBATIM
  return 0;
  ENDVERBATIM
}

FUNCTION exptable(x) { 
  TABLE  FROM -10 TO 10 WITH 2000
  
  if ((x > -10) && (x < 10)) {
    exptable = exp(x)
  } else {
    exptable = 0.
  }
}

PROCEDURE rates(v(mV)) {
  TABLE B
  DEPEND mg
  FROM vmin TO vmax WITH 200
  : Stevens & Jahr 1990a,b
  
  B = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}
