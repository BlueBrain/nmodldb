TITLE passive membrane channel

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
}

NEURON {
	SUFFIX pasR
	NONSPECIFIC_CURRENT i
	RANGE g, e, vol
}

PARAMETER {
	g = .001	(S/cm2)	<0,1e9>
	e = -70	(mV)
}

ASSIGNED {v (mV)  i (mA/cm2)  vol (mV)}

INITIAL {
  vol = v-e
}

BREAKPOINT {
	i = g*(v - e)
        vol=v-e
}
