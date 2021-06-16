NEURON {
	POINT_PROCESS gap
	NONSPECIFIC_CURRENT i
	RANGE r, i
	POINTER vgap
}
PARAMETER {
		
	r = 2000 (megohm)
}
ASSIGNED {
	i (nanoamp)
	vgap (millivolt)
	v (millivolt)
}
BREAKPOINT {
	i = (v - vgap)/r
}

