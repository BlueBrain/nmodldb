NEURON {POINT_PROCESS AMPA}

PARAMETER {
  Cdur	= 1.1	(ms)		: transmitter duration (rising phase)
  Alpha	= 10	(/ms mM)	: forward (binding) rate
  Beta	= 0.5	(/ms)		: backward (unbinding) rate
  Erev	= 0	(mV)		: reversal potential
}
INCLUDE "netcon.inc"
:** NMDA
