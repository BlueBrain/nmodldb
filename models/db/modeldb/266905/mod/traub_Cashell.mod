TITLE Calcium dynamics Traub

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
    (mS) = (millisiemens)
   	(molar) = (1/liter)
	(uM)    = (micromolar)
}

NEURON {
    THREADSAFE
    : note - every variable accessible in NEURON will be having the suffix _CaShell
    
    SUFFIX CaShell
    USEION ca READ ica WRITE cai
    RANGE phi
    GLOBAL beta_x, cai0
}

PARAMETER {
    phi = 578.135  (uM-cm2/mA-ms) : value at soma
    beta_x = 0.075	(/ms)
    cai0 = 0	(mM)
}

STATE {
    cai	(uM)
}

ASSIGNED {
    ica	(mA/cm2)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
}

INITIAL {
   cai = cai0
}

DERIVATIVE states {  
    cai' = - phi*ica -beta_x*cai
}

