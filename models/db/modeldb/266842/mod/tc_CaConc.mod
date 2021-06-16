TITLE Thalamocortical neuron calcium concentration

COMMENT

	Calcium concentration
	Implementation of Meijer et al., 2011
	Written by Xu Zhang, UConn, 2018

ENDCOMMENT

NEURON {
    SUFFIX tcCaConc
    USEION ca READ ica WRITE cai
    RANGE cai, kCa, dep
}

UNITS {
    (molar) = (1 / liter)
    (mM) = (millimolar)
    (mV) = (millivolt)
    (mA) = (milliamp)
}

PARAMETER {
	cabuf	= 2.4e-4 (mM)
    tau_ca = 5.0 (ms)
	kca = 5.1821e-5 (1/coulomb)
    depth = 0.2 (micron)
}

ASSIGNED {
    ica (mA/cm2)
}

STATE {
    cai (mM)
}

INITIAL {
    cai = cabuf
}

BREAKPOINT {
    SOLVE states METHOD cnexp
}

DERIVATIVE states {
    cai' = (cabuf-cai)/tau_ca - kca*ica
}
