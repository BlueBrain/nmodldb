TITLE Leak Current

COMMENT

NEURON implementation of a passive leak/shunt

Laboratory for Neuronal Circuit Dynamics
RIKEN Brain Science Institute, Wako City, Japan
http://www.neurodynamics.brain.riken.jp

Date of Implementation: April 2005
Contact: akemann@brain.riken.jp

ENDCOMMENT

NEURON {
	SUFFIX leak
	NONSPECIFIC_CURRENT i
	RANGE i, e, gbar
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(nA) = (nanoamp)
	(pA) = (picoamp)
	(S)  = (siemens)
}

PARAMETER {
	gbar = 5e-5 (S/cm2)  < 0, 1e9 > : 5e-5 mho/cm2 = 20.000 ohm cm2
	e = -60 (mV)
}

ASSIGNED {
	i (mA/cm2)
	v (mV)
}

BREAKPOINT {
	i = gbar*(v - e)
}
