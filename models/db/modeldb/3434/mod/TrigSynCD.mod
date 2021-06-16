TITLE TrigSynCD.mod   Triggerable alpha synapse

COMMENT
 alpha function synapse
 written by Jonathan Z Simon, jzsimon@isr.umd.edu
 but core code taken from Hines & Carnivale, "Expanding NEURON with NMODL"
ENDCOMMENT

COMMENT
 modified by Jonathan Simon, jzsimon@isr.umd.edu
 a synaptic current with alpha function conductance that works out to be
  i = g * (v - e)      i(nanoamps), g(micromhos);
 where
  g = 0 before the synapse is triggered at onset
  g = stim * (t - onset)/tau * exp(-(t - onset - tau)/tau) for t > onset
 this has the property that the maximum value is gmax and occurs at
  t =  onset + tau.
 Except that this implementation is more flexible and allows multple
  overlapping inputs and corresponding summing conductances automatically.
ENDCOMMENT
					       
NEURON {
	POINT_PROCESS TrigSynCD
	RANGE tau, gmax, e, i, g
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau=.1 (ms)
	gmax=0 	(uS)
	e=0	(mV)
}

ASSIGNED {
	i (nA)
	bath (uS)
	g (uS)
	k (/ms)
	v	(mV)
}

CONSTANT { exp1 = 2.7182818284590452354} :Euler constant for correct max peak

STATE { A (uS) G (uS) }

INITIAL {
	k = 1/tau
	A = 0
	G = 0
}


BREAKPOINT {
	SOLVE state METHOD sparse
	i = G*(v - e)
	g = G
}

KINETIC state {
	~ A <-> G	(k, 0)
	~ G -> 	(k)
}

NET_RECEIVE (weight (uS)) {
    state_discontinuity(A, A + gmax*exp1)
}
