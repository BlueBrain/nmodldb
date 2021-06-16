NEURON {
	SUFFIX leak_gp
	NONSPECIFIC_CURRENT i
	RANGE i, e, g, gbar
	GLOBAL noisevar
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
}

PARAMETER {
	gbar = 1	(S/cm2)
	e = 0		(mV)
	noisevar=0
}

ASSIGNED {
	g	(S/cm2)
	i 	(mA/cm2)
	v	(mV)
}

BREAKPOINT {
	g = gbar*normrand(1,noisevar)
	i = g*(v-e)
}