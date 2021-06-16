INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
NEURON {
	POINT_PROCESS gapvdep
	NONSPECIFIC_CURRENT i
	RANGE gg,i,slow
	POINTER vpregap
}
PARAMETER {
	v (millivolt)
	vpregap (millivolt)
	gg = 5e-4 (micromho)
}

STATE { slow
}
ASSIGNED {
	i (nanoamp)
}
BREAKPOINT {
	SOLVE sssss METHOD cnexp
	i = (v - vpregap)*gg*(vdepfac(v-slow))
VERBATIM

ENDVERBATIM
}

DERIVATIVE sssss {
	slow' = (vpregap-slow)/tau_slow(vpregap)
}




FUNCTION vdepfac(v){
vdepfac = 1.0/(1+exp(-(0.12*(v-5))))

}

FUNCTION tau_slow(ww){
tau_slow=0.2+200/(1+exp(-0.5*(ww+50)))
}







