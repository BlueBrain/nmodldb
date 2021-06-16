TITLE Ca influx based on ica, cai

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (uA) = (microamp)
    (molar) = (1/liter)
    (mM) = (millimolar)
    (uM) = (micromolar)
}

NEURON {
	THREADSAFE
    SUFFIX CaInternal
    USEION ca READ ica, cai WRITE cai
    GLOBAL ca_min
    RANGE alpha_ca, tau_ca, ca_init
}

PARAMETER {
: parameters can be set in hoc template files
    :alpha_ca = 0.006 (mM-cm2/ms/mA)
    alpha_ca = 0.003 (uM-cm2/ms/uA) 
    tau_ca = 500 (ms)
    ca_init=5e-5 (mM)
    ca_min=1e-5	(mM)
}

ASSIGNED {
    ica (mA/cm2)
}

STATE {
    cai (mM)
}

BREAKPOINT {
    SOLVE state METHOD derivimplicit
}

DERIVATIVE state {
    cai' = -1*alpha_ca*ica - ((cai-ca_min)/tau_ca)
}

INITIAL {
    cai = ca_init
}






