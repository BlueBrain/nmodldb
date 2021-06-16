NEURON {
        POINT_PROCESS GrCtonicGABA
        RANGE  e, i
        NONSPECIFIC_CURRENT i

        RANGE g
}

UNITS {
       	(nA) = (nanoamp)
        (mV) = (millivolt)
        (uS) = (microsiemens)
}

PARAMETER {
	        e=-75   (mV)
	        g=260e-6 (uS)

}

ASSIGNED {
    v (mV)
	i (nA)

}


INITIAL {


            i = g*(v - e)
}

BREAKPOINT {

    i = g*(v - e)

}
