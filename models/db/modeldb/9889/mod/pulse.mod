: $Id: pulse.mod,v 1.2 1998/07/21 17:23:37 billl Exp $
COMMENT
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  POINT_PROCESS PULSE
  NONSPECIFIC_CURRENT i
  RANGE amp, dur
}
 
UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
  (umho) = (micromho)
  (mM) = (milli/liter)
}

PARAMETER {
  DELAY = 0                     : axonal delay ??
  dur	= 0	(ms)            : duration
  amp   = 0     (nA)            : amplitude
}

ASSIGNED {
  Aon
  i 		(nA)
  dt
}

INITIAL {
  i = 0
}

BREAKPOINT {
  i = Aon
}


NET_RECEIVE(weight, on, nspike) {
  if (flag == 0 && !on) { : turn on
    nspike = nspike + 1 : why?
    Aon = -amp
    : ?? state_discontinuity() : need to signal an at_time but no states to reset
    net_send(dur, nspike)
  }
  if (flag == nspike) { : if this associated with last pulse then turn off
    Aon = 0
    on = 0
  }
}
