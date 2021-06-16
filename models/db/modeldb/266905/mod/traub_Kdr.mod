TITLE delayed rectifier K channel Traub

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
    (mS) = (millisiemens)
}

NEURON {
    THREADSAFE
    : note - every variable accessible in NEURON will be having the suffix _gKdr
    
    SUFFIX gKdr
    USEION k READ ek WRITE ik
    RANGE gmax, g, i, tau_n, n_inf
    : GLOBAL n_inf, tau_n
}

PARAMETER {
    gmax=0.015 (mho/cm2)  : value at soma
}

STATE {
    n
}

ASSIGNED {
    v (mV)
    ek (mV)
    ik		(mA/cm2)
    n_inf	(1)
    tau_n	(ms)
    g		(S/cm2)
    i		(mA/cm2)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g = gmax*n
    i = g*(v-ek)
    ik = i
}

INITIAL {
    rates(v)
    n = n_inf
}


FUNCTION alpha_n(v(mV)) {
    alpha_n = 0.016*(35.1-v)/(exp( (35.1-v)/5 ) - 1)
}

FUNCTION beta_n(v(mV)) {
    beta_n = 0.25*exp( (20 - v)/40 )
}


DERIVATIVE states {  
    rates(v)
    n' = (n_inf - n)/tau_n
}

PROCEDURE rates(v (mV)) { :callable from hoc
    LOCAL alpha, beta
    TABLE n_inf, tau_n  : DEPEND vhalfn, vhalfl, tlmax, tnmax
    FROM -20 TO 130 WITH 750
    
    alpha = alpha_n(v)
    beta = beta_n(v)
    tau_n = 1/(alpha + beta)     
    n_inf = alpha*tau_n
    
}

