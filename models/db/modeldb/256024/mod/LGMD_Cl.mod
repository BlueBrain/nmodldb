TITLE Cl influx based on icl, cli

UNITS {
    (mV) = (millivolt)
	(mA) = (milliamp)
    (uA) = (microamp)
    (molar) = (1/liter)
    (mM) = (millimolar)
}

NEURON {
    SUFFIX ClInternal
    USEION cl READ icl, cli WRITE cli VALENCE -1
    RANGE alpha_cl, tau_cl
}

PARAMETER {
	alpha_cl = 0.006 (mM-cm2/ms/mA)          : both parameters can be set in hoc template files
    tau_cl = 50 (ms)
    cl_init = 5 (mM)
    clo = 157 (mM)
}

ASSIGNED {
    icl (mA/cm2)
}

STATE {
    cli (mM)
}

BREAKPOINT {
    SOLVE state METHOD derivimplicit
}

DERIVATIVE state {
    cli' = -1*alpha_cl*icl - (cli/tau_cl)
}

INITIAL {
    cli = cl_init
}












