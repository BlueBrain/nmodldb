TITLE passive membrane channel

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
}

NEURON {
	SUFFIX pas_gp
	NONSPECIFIC_CURRENT i
	RANGE g, e, gbar
}

PARAMETER {
	gbar = 1	(S/cm2)
	e = -70		(mV)
}

ASSIGNED {
	v	(mV)
	i	(mA/cm2)
	g	(S/cm2)
}

BREAKPOINT {
	g = gbar
	i = g*(v - e)
}
