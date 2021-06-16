NEURON {
	POINT_PROCESS Gap
	POINTER vnb
	RANGE r, i
	NONSPECIFIC_CURRENT i
}
PARAMETER { r = 1e10 (megohm) }
ASSIGNED {
	v 		(millivolt)
	vnb 	(millivolt)
	i 		(nanoamp)
}
BREAKPOINT { i = (v - vnb)/r }