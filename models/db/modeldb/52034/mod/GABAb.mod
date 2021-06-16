NEURON {  POINT_PROCESS GABAb }

PARAMETER {
  Cdur	= 85	(ms)		: transmitter duration (rising phase)
  Alpha	= 0.016	(/ms mM)	: forward (binding) rate
  Beta	= 0.0047 (/ms)		: backward (unbinding) rate
  Erev	= -90	(mV)		: reversal potential
  DELAY = 0
  Deadtime = 1. (ms)		: mimimum time between release events
  GMAX	= 1.0	(uS)		: maximum conductance
  Thresh = -70
}

INCLUDE "sns.inc"

    : EXTRA BREAKPOINT MUST BE BELOW THE INCLUDE
    BREAKPOINT {
      if (v>Thresh) {
        g = g*((v/Thresh)*(v/Thresh)*(v/Thresh))
        i = g*(v - Erev)
      }
    }
