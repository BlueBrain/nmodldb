INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
NEURON {
	POINT_PROCESS gap
	NONSPECIFIC_CURRENT i
	RANGE gg, i
	POINTER vgap
}
PARAMETER {
	v (millivolt)
	vgap (millivolt)
	gg = 5e-4 (micromho)
}
ASSIGNED {
	i (nanoamp)
}
BREAKPOINT {
	i = (v - vgap)*gg
}
