TITLE Leak for horizontal cell
: Leak current for horizontal cells
: 
: Based on parameters of Aoyama et al. (2000)


NEURON {
	SUFFIX HzLeak
    NONSPECIFIC_CURRENT i
    RANGE i, e, gbar
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
    (mS) = (millimho)
}

PARAMETER {
    gbar = 4.7170e-06 (mho/cm2) : 0.5 nS total
    e = -80 (mV)
}

ASSIGNED {
	i    (mA/cm2)
    v    (mV)
}

BREAKPOINT {
    i = gbar * (v - e)
}

