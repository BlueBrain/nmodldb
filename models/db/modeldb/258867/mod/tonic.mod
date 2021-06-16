TITLE tonic GABA current according to Pavlov et al., 2009 in J Neurosci


NEURON{
    SUFFIX tonic
    NONSPECIFIC_CURRENT itonic
    RANGE  itonic, e, gtonic, minf
}


UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
    gtonic=0.001 (siemens/cm2) <0, 1e9>
    e = -70 (mV)
}

ASSIGNED{
    v (mV)
    itonic (milliamp/cm2)
    minf     
}

BREAKPOINT {
    minf=a(v)/(a(v)+b(v))
    itonic=gtonic*minf*(v-e)
}

FUNCTION a(v(mV)) {
    LOCAL x
    if (fabs(v+20) > 1e-5) {
        x = 0.1 * (v + 20)
    }else{
        x = 0.1
    }
    a = (50 * x / (1 - exp(-x)))
}

FUNCTION b(v(mV)) {
    LOCAL x
    if (fabs(v-10) > 1e-5) {
        x = -0.08 * (v - 10)
    } else {
        x = -0.08
    }
    b = (20 * x / (1 - exp(-x)))
}