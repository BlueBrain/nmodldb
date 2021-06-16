: $Id: pulse.mod,v 1.1 1998/07/21 13:54:10 billl Exp $
COMMENT
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  POINT_PROCESS PULSE
  NONSPECIFIC_CURRENT i
  RANGE amp, dur
}
 
INCLUDE "snsarr.inc"  : array management routines

UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
  (umho) = (micromho)
  (mM) = (milli/liter)
}

PARAMETER {
  dur	= 0	(ms)            : duration
  amp   = 0     (nA)            : amplitude
  Cdur  = 0                     : not used
  Deadtime = 0	(ms)		: mimimum time between release events, not used
  GMAX = 1	(umho)		: not used
  DELAY = 0     (ms)
}

ASSIGNED {
  i 		(nA)
  dt
}

INITIAL {
  i = 0
  if (nsyn > 0) {
    initq()   
  } 
}

BREAKPOINT {
  if (nsyn>0) { : do not try accessing a queue that hasn't been allocated
  VERBATIM 
  static int ii,who;
  static QueU *pqueu;
  static SynS *ppst;

  pqueu = (QueU *)((unsigned long) queu);

  while (t >= pqueu[(int)begsyn].time) { /*  somebody spiked delay time ago */
    i = -amp;
    popqh1(dur);		/* next (also add Cdur to value on the queu) */
  }

  while (t >= pqueu[(int)endsyn].time) { /*  somebody needs to be turned off */
    i = 0;
    popqh2();  /* next */
  }
  ENDVERBATIM
}
}
