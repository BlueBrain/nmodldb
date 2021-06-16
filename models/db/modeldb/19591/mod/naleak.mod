TITLE Ca leak

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

INDEPENDENT { v FROM -100 TO 50 WITH 50	(mV) }

NEURON {
	SUFFIX naleak
	NONSPECIFIC_CURRENT i
	RANGE gbar, i
        GLOBAL e
}

PARAMETER {
	gbar = 1.167e-5	(mho/cm2)
	e = 55	(mV)
}

ASSIGNED { i    (mA/cm2) }

BREAKPOINT {
	i = gbar*(v - e)
}
