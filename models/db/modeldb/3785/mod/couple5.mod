
TITLE couple5.mod

NEURON {
	SUFFIX couple5
	NONSPECIFIC_CURRENT i
	RANGE r, i
	POINTER vc
}

PARAMETER {
	r = 100 (ohm-cm2)
}

ASSIGNED {
	i (milliamp/cm2)
	v (millivolt)
	vc (millivolt)
}

BREAKPOINT {
	i = (v - vc)/r
}
