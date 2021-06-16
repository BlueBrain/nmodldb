NEURON {
	POINT_PROCESS gap
	NONSPECIFIC_CURRENT i
	RANGE del, r, i
	POINTER vgap
}
PARAMETER {
	v (millivolt)
	vgap (millivolt)
	r = 1e10 (megohm)
	del(ms)
}
ASSIGNED {
	i (nanoamp)
}
BREAKPOINT {
	if (t > del) {
		i = (v - vgap)/r
	}else{
		i=0
	}
}
