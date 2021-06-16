NEURON {POINT_PROCESS AMPA}

PARAMETER {
  Cdur	= 1.1	(ms)		: transmitter duration (rising phase)
  Alpha	= 1	(/ms mM)	: forward (binding) rate
  Beta	= 0.5	(/ms)		: backward (unbinding) rate
  Erev	= 0	(mV)		: reversal potential
  Deadtime = 2.5 (ms)		: mimimum time between release events
  GMAX     = 1	(uS)		: maximum conductance
  DELAY = 0
}
INCLUDE "netcon.inc"
:* DEFAULTS
