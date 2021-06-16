: constant nmda conductance
: from Durstewitz & Gabriel (2006), Cerebral Cortex

NEURON {
	SUFFIX nmdac
	NONSPECIFIC_CURRENT inmdac
        RANGE gNMDAcbar
}

UNITS {
        (uS) = (microsiemens)
        (nA) = (nanoamp)
        (mV) = (millivolt)
}

PARAMETER {
	gNMDAcbar=0.0 (mho/cm2) <0,1e9>
	enmda=0.0     (mV)
}

ASSIGNED {
	v (mV)
	inmdac
}

BREAKPOINT {
	LOCAL s
	s = 1.50265/(1+0.33*exp(-0.0625*v))
	inmdac = gNMDAcbar*s*(v-enmda)
}
