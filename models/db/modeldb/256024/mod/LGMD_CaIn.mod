TITLE Ca influx based on ica, cai
: adaptation of LGMD_Ca to include eflux 2 sources
: 1 low volume high affinity, the other high volume low afinity

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
    SUFFIX CaIn
    USEION ca READ ica, cai WRITE cai
    GLOBAL ca_min1, ca_min2, ca_init
    RANGE tau1, tau2, alpha_ca
}

PARAMETER {
: parameters can be set in hoc template files
    :alpha_ca = 0.006 (mM-cm2/ms/mA)	: sets influx
    alpha_ca = 0.001 (uM-cm2/ms/uA) 
    tau1 = 400		(ms)
    tau2 = 20		(ms)
    ca_init = 1e-4	(mM)
    ca_min1 = 2e-5	(mM)
    ca_min2 =1.4e-4	(mM)
}

ASSIGNED {
    ica (mA/cm2)
    eflux (mM/ms)
}

STATE {
    cai (mM)
}

BREAKPOINT {
    SOLVE state METHOD derivimplicit
}

DERIVATIVE state {
	if (cai>ca_min2) {
	   eflux = (cai-ca_min1)/tau1 + (cai-ca_min2)/tau2
	} else if (cai>ca_min1) {
	    eflux = (cai-ca_min1)/tau1
	} else {
		eflux = 0
	}
	cai' = -1*alpha_ca*ica - eflux
}

INITIAL {
    cai = ca_init
}





