: leak.mod codes a leak (passive) channel.
:
: Takaki Watanabe
: wtakaki@m.u-tokyo.ac.jp

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

NEURON {
	SUFFIX leak
	NONSPECIFIC_CURRENT i
	RANGE g, erev
}

PARAMETER {
	g = .001	(mho/cm2)
	erev = -85	(mV)
}

ASSIGNED {
i(mA/cm2)
v(mV)
}

BREAKPOINT {
	i = g*(v - erev)
}


