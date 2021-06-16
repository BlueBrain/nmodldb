NEURON{ POINT_PROCESS NMDA
  RANGE B 
}

PARAMETER {
  mg    = 1.    (mM)            : external magnesium concentration
  Cdur	= 1.	(ms)		: transmitter duration (rising phase)
  Alpha	= 4.	(/ms mM)	: forward (binding) rate
  Beta	= 0.0067 (/ms)		: backward (unbinding) rate 1/150
  Erev	= 0.	(mV)		: reversal potential
  Deadtime = 1	(ms)		: mimimum time between release events
  GMAX	= 1     (S/cm2)         : maximum conductance
  DELAY = 0                     : axonal delay
}

ASSIGNED { B }
INCLUDE "sns.inc"
: BREAKPOINT MUST GO AFTER INCLUDE
BREAKPOINT {
  rates(v)
  g = g * B
  i = i * B
}
PROCEDURE rates(v(mV)) {
  TABLE B
  DEPEND mg
  FROM -100 TO 80 WITH 180
  B = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}

:* >>>> Utility defaults <<<<  
