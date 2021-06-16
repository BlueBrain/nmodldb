INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS int1_lgsyn
	POINTER vpre
	RANGE gmax, g, e, i
	NONSPECIFIC_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	gmax=0 	(umho)
	e=0	(mV)
	v	(mV)
}

STATE { synon synoff
}

ASSIGNED {
	i (nA)
	g (umho)
	vpre (mV)
}



BREAKPOINT {

	SOLVE syn METHOD sparse
	g = gmax *synon
	i = g*(v - e)

}

KINETIC syn {

	~ synoff <-> synon (syninf(vpre)/tausyn(vpre),(1-syninf(vpre))/tausyn(vpre))

}

INITIAL {
	synon = 0.0
	synoff = 1.0
}

FUNCTION syninf(v){
	syninf = 1/(1+exp(-0.5*(v+49)))
}

FUNCTION tausyn(v){
:	tausyn = 2+98/(1+exp(-0.5*(v+49)))
	tausyn = 50
}
