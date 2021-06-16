NEURON {  POINT_PROCESS GABAA }
PARAMETER {
  Cdur	= 1.0	(ms)		: transmitter duration (rising phase)
  Alpha	= 0.53	(/ms mM)	: forward (binding) rate
  Beta	= 0.18	(/ms)		: backward (unbinding) rate
  Erev	= -80	(mV)		: reversal potential
}
INCLUDE "netcon.inc"
:** GABAB
