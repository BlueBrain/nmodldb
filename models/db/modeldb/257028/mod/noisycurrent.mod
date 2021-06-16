TITLE Membrane noise

NEURON {
	POINT_PROCESS NoisyCurrent
	RANGE i, noise
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
}

PARAMETER {
	noise = 0 (nA)
}

ASSIGNED {
	i (nA)
}

BREAKPOINT {
	i = noise
}