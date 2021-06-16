TITLE Ca influx based on ica, cai

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX CaInternal
    USEION ca READ ica, cai WRITE cai
    RANGE alpha_ca, tau_ca
}

PARAMETER {
    alpha_ca = .006 (mM*cm2/ms/mA)          : both parameters are set in hoc template files
    tau_ca = 130 (ms)
}

ASSIGNED {
    ica (mA/cm2)
}

STATE {
    cai
}

BREAKPOINT {
    SOLVE state METHOD derivimplicit
}

DERIVATIVE state {
    cai' = -1*alpha_ca*ica - (cai/tau_ca)
}

INITIAL {
    cai = 0
}
