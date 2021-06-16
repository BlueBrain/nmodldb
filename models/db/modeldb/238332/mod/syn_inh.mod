NEURON {
    POINT_PROCESS syn_inh
    RANGE tau_o                           : parameter
    RANGE tau_c                          : parameter
    RANGE erev                              : parameter
    RANGE i                                 : exposure
    RANGE syn_step
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
    tau_o = 1.5 (ms)
    tau_c = 4.0 (ms)
    erev = -75 (mV)
    syn_step = 3.0
}

ASSIGNED {
    v (mV)
    i (nA)
}

STATE {
    o (1/ms)
    c (1/ms)
}

INITIAL {
    o = 0
    c = 0
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    i = (c - o) * (v-erev)
}

NET_RECEIVE(weight (uS)) {
    :printf("w = %g \n", syn_step*weight)
    o = o + syn_step*weight
    c = c + syn_step*weight
}

DERIVATIVE states {
    o' = -o/tau_o
    c' = -c/tau_c
}
