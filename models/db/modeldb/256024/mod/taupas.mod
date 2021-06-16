TITLE voltage independent (passive) channel with slowly changing conductance (introduces inductance)

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S) = (siemens)
}

NEURON {
    THREADSAFE
    SUFFIX taupas

    NONSPECIFIC_CURRENT i
	RANGE g,tau, i
}

PARAMETER {
    g = 0.0001	(S/cm2)	<1e-8,10>
    e = -65		(mV)
	tau = 30	(ms)	<0,1e4>
}

ASSIGNED { 
    v	(mV)
    i	(mA/cm2)
}

STATE {
	lv	(mV)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    if (tau<1e-3) {
		i = g*(v-e)
    } else {
		i = g*(lv-e)
    }
}

INITIAL {
    lv = v
}

DERIVATIVE states {
	lv' = (v - lv)/tau
}



