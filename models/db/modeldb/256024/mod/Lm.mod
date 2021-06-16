TITLE membrane inductance

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S) = (siemens)
}

NEURON {
    THREADSAFE
    SUFFIX Lm

    NONSPECIFIC_CURRENT i

	RANGE L,maxg, i
}

PARAMETER {
	maxg = 4e-5	(S/cm2)
    e = -65		(mV)
	L = 2.0e3	(henry-cm2)
}

ASSIGNED { 
	dt	(ms)
    v	(mV)
    i	(mA/cm2)
}

STATE {
	lv	(mV)
}

BREAKPOINT {
	:tauL = L*maxg*1e3(ms/s)
	if (L>0) {
		i = i + (v-e)*dt*1e-3(s/ms)/L
	}
	if ( (v>e && i>maxg*(v-e)) ||  (v<e && i<maxg*(v-e)) ) { i=maxg*(v-e) }
}

INITIAL {
    i = 0
}



