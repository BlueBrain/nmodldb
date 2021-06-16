TITLE passive membrane channel

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
}

NEURON {
	THREADSAFE
	SUFFIX pasnts
	NONSPECIFIC_CURRENT i
	RANGE g, e

}

PARAMETER {
	g = .00001	(S/cm2)	<0,1e9>
	e = -50	(mV)
}

ASSIGNED {v (mV)  i (mA/cm2)}

BREAKPOINT {
	i = g*(v - e)
}
