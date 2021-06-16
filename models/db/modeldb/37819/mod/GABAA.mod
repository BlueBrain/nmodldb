NEURON {  POINT_PROCESS GABAA }
PARAMETER {
  Cdur	= 0.3	(ms)		: transmitter duration (rising phase)
  Alpha	= 10	(/ms mM)	: forward (binding) rate
  Beta	= 0.16	(/ms)		: backward (unbinding) rate
  Erev	= -80	(mV)		: reversal potential
}
INCLUDE "netcon.inc"
:** GABAB2
