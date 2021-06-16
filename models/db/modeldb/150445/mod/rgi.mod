TITLE Random current


NEURON {
    SUFFIX rgi
    RANGE dc, sd, driver
    NONSPECIFIC_CURRENT i
}

UNITS {
    (mA) = (milliamp)
    (mA/cm2) = (nanoamp/cm2)
}


PARAMETER {
	dt	     (ms)
	dc	= 0. (mA/cm2)	: DC offset of the overall current
	sd	= 0. (mA/cm2)	: square root of the steady-state variance of the current
}

ASSIGNED {
    i (mA/cm2)              : overall sinusoidal noisy current
    driver
}

INITIAL {
    i = dc
}

BREAKPOINT {
    i = dc - sd*driver
}





