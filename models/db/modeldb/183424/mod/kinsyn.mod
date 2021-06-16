COMMENT
Kinetic model of synaptic transmission
dr/dt = (rinf(T)-r)/taur(T)

i = g * r * ( v - e )
ENDCOMMENT

NEURON {
    POINT_PROCESS KineticSynapse
    RANGE onset, duration, tau1, tau2, gmax, e, i
    NONSPECIFIC_CURRENT i
}
UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (uS) = (microsiemens)
}

PARAMETER {
    onset=0 (ms)
    duration=1 (ms)
    tau1=0.1 (ms)
    tau2=5 (ms)
    gmax=0 (uS)
    e=0 (mV)
}

ASSIGNED {
    v (mV)
    i (nA)
    g (uS)
}

STATE {
    r
}

BREAKPOINT {
    if (gmax) {
	at_time( onset )
	at_time( onset + duration )
    }
    SOLVE state METHOD cnexp
    i = gmax * r * ( v - e )
}

UNITSOFF

INITIAL {
    r = 0
}

DERIVATIVE state {
    r' = ( 1.0 / tau1 - 1.0 / tau2 ) * T( t ) * ( 1.0 - r ) - r / tau2
}

FUNCTION T( t ) {
    if ( t < onset || t > onset + duration ) {
	T = 0
    } else {
	T = 1
    }
}

UNITSON
