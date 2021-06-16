TITLE Leak Density Mechanism

NEURON {
	SUFFIX leak
	NONSPECIFIC_CURRENT i
	RANGE i, e, gmax
}

PARAMETER {
	gmax = 0.00013333 (siemens/cm2) 
	e = -65 (millivolt)
}

ASSIGNED {
	i (milliamp/cm2)
	v (millivolt)
}

BREAKPOINT { 
	i = gmax * (v - e) 
}
