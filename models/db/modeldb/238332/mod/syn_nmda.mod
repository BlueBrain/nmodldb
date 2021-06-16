NEURON {
    POINT_PROCESS syn_nmda
    RANGE tau_o                           : parameter
    RANGE tau_c                          : parameter
    RANGE erev                              : parameter
    RANGE c1
    RANGE c2
    RANGE syn_step
    RANGE i                                 : exposure
    NONSPECIFIC_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
    (uA) = (microamp)
    (mA) = (milliamp)
    (A) = (amp)
    (mV) = (millivolt)
    (mS) = (millisiemens)
    (uS) = (microsiemens)
    (molar) = (1/liter)
    (kHz) = (kilohertz)
    (mM) = (millimolar)
    (um) = (micrometer)
    (S) = (siemens)
}

PARAMETER {
    tau_o = 5.0 (ms)
    tau_c = 80.0 (ms)
    erev = 0.0 (mV)
    c1 = 0.05
    c2 = -0.08
    syn_step = 1.25
}

ASSIGNED {
    v (mV)
    i (nA)
}

STATE {
    o
    c
}

INITIAL {
    o = 0
    c = 0
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    i = mgBlock(v) * (c - o) * (v-erev)
}

NET_RECEIVE(weight (uS)) {
    o = o + syn_step*weight
    c = c + syn_step*weight
}


DERIVATIVE states {
    o' = -o/tau_o
    c' = -c/tau_c
}

UNITSOFF
FUNCTION mgBlock(v (mV)) {
    mgBlock = 1 / (1 + c1*exp(c2*v))
}
UNITSON
