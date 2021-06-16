: A synthetic conductance

NEURON {
	POINT_PROCESS Gsynth
	NONSPECIFIC_CURRENT i
	RANGE i, e, g
}

PARAMETER {
	g = 0 (microsiemens) < 0, 1e9 > : g >= 0
	e = 0 (millivolt)
}

ASSIGNED {
	i (nanoamp)
	v (millivolt)
}

BREAKPOINT { i = g*(v-e) }
