TITLE Lpas channel for LGMD

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S) = (siemens)
}

NEURON {
    THREADSAFE
    SUFFIX Lpas

    NONSPECIFIC_CURRENT i

	RANGE g0, g,tauL,pl,L,e, i
}

PARAMETER {
    g0 = 0.0001	(S/cm2)
    e = -65		(mV)
	pl = 0.4	(1)	<0,1>
	L = 1.0e3	(henry-cm2)
}

ASSIGNED { 
    v	(mV)
    i	(mA/cm2)
    tauL(ms)
    g	(S/cm2)
}

STATE {
	lv	(mV)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
	tauL = L*g0*pl*1e3(ms/s)
    if (tauL<1e-2)	{ tauL = 1e-2 (ms) }
:	i = g0*(v-e) + g0*pl*(lv-e)
	i = g0*(v-e+pl*(lv-v))
    if (fabs((v-e)/1(mV))>1e-3) {
		g = i/(v-e)
	}
}

INITIAL {
    lv = v
    g=0
}

DERIVATIVE states {
	: if L is very low and pl is not, there are problems with high freq oscillations
	lv' = (v - lv)/tauL
}



