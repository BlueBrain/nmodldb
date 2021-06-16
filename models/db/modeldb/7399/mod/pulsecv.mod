: $Id: pulsecv.mod,v 1.9 1998/10/16 20:33:58 billl Exp $
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
  dur	= 0	(ms)            : duration
  amp   = 0     (nA)            : amplitude
}

ASSIGNED {
  on
  i 		(nA)
  dt
}

INITIAL {
  on = 0
  i = 0
}

BREAKPOINT {
  if (on==1) { 
    i = -amp
  } else {
    i = 0
  }
}

NET_RECEIVE(weight) {   : flag is an implicit argument and  normally 0
  if (flag == 0) { : a spike, so turn on 
    on = 1
    net_send(dur,-1) : come again in dur with flag set
  }
  if (flag == -1) { : turn off
    on = 0
  }
}
