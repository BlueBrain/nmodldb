NEURON {  POINT_PROCESS AMPA }

PARAMETER {
  Cdur	= 1	(ms)		: transmitter duration (rising phase)
  Alpha	= 1.	(/ms mM)	: forward (binding) rate
  Beta	= 0.5	(/ms)		: backward (unbinding) rate
  Erev	= 0	(mV)		: reversal potential
  DELAY = 0
  Deadtime = 1. (ms)		: mimimum time between release events
  GMAX	= 1.0	(uS)		: maximum conductance
}

INCLUDE "sns.inc"
