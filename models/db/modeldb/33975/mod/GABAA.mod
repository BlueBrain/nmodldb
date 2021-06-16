NEURON {  POINT_PROCESS GABAA }
PARAMETER {
  Cdur	= 0.3	(ms)		
  Alpha	= 12	(/ms mM)	
  Beta	= 0.1	(/ms)		
  Erev	= -75	(mV)		
}
INCLUDE "netcon.inc"

:** GABAB2
