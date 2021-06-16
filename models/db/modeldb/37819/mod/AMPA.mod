NEURON {POINT_PROCESS AMPA}

PARAMETER {
  Cdur	= 0.3	(ms)		: transmitter duration (rising phase)
  Alpha	= 0.47	(/ms mM)	: forward (binding) rate
  Beta	= 0.18	(/ms)		: backward (unbinding) rate
  Erev	= 0	(mV)		: reversal potential
}
INCLUDE "netcon.inc"
:** NMDA
