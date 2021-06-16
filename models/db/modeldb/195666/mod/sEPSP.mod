TITLE sEPSP

COMMENT
-----------------------------------------------------------------------------

    sEPSP current model for summation analysis
    ==================================================

 IMPLEMENTATION

  This mechanism is implemented as a nonspecific current defined as a
  point process, mimicking a current-clamp stimulation protocol, injecting
  simulated EPSP I(t).

  I = 0 for t < onset and
  I(t) = A*(1 - exp(-1*(t-onset)/tau_r)) * (exp(1 - (t-onset) / tau_f))
	  for t > onset > t+offset
  I = 0 for t > offset


 PARAMETERS

  This mechanism takes the following parameters:

  A		= 1.0 (nA) : initial current offset
  taur	= 0.3 (ms) : time constant of rise
  tauf	= 3.0 (ms) : time constant of fall
  onset	= 0.0 (ms) : start time of current
  offset = onset+20 (ms)
-----------------------------------------------------------------------------
ENDCOMMENT


NEURON {
    POINT_PROCESS sEPSP
    RANGE A, onset, i
    NONSPECIFIC_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
}

PARAMETER {
	A		= 0.0 (nA)	: initial current offset.
	taur	= 0.3 (ms)	: time constant of rise.
	tauf	= 3.0 (ms)	: time constant of fall.
	onset	= 0.0 (ms)	: start time of current
}

ASSIGNED {
    i     (nA)        : injected current
}


BREAKPOINT {

    if ((t < onset) || (t > onset+20)) {
    	i = 0
	} else { 
		i = -A*(1 - exp(-1*(t-onset)/taur)) * (exp(1 - (t-onset) / tauf))
    }
}


