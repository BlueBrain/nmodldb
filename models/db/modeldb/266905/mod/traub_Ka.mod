TITLE A-type K channel Traub

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
    (mS) = (millisiemens)
}

NEURON {
    THREADSAFE
    : note - every variable accessible in NEURON will be having the suffix _gKa
    
    SUFFIX gKa
    USEION k READ ek WRITE ik
    RANGE gmax, g, i
	: GLOBAL a_inf, tau_a, b_inf, tau_b
}

PARAMETER {
    gmax=0.005 (mho/cm2)  : value at soma
}

STATE {
    a
    b
}

ASSIGNED {
    v (mV)
    ek (mV)
    ik		(mA/cm2)
    a_inf	(1)
    b_inf	(1)
    tau_a	(ms)
    tau_b	(ms)
    g		(S/cm2)
    i		(mA/cm2)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g = gmax*a*b
    i = g*(v-ek)
    ik = i
}

INITIAL {
    rates(v)
    a = a_inf
    b = b_inf
}


FUNCTION alpha_a(v(mV)) {
    alpha_a = 0.02*(13.1-v)/(exp( (13.1-v)/10 ) - 1)
}

FUNCTION beta_a(v(mV)) {
    beta_a = 0.0175*(v - 40.1)/(exp( (v-40.1)/10 ) - 1)
}

FUNCTION alpha_b(v(mV)) {
    alpha_b = 0.0016*exp( (-13-v)/18 )
}

FUNCTION beta_b(v(mV)) {
    beta_b = 0.05/(1 + exp( (10.1 - v)/5 ) ) 
}

DERIVATIVE states {  
    rates(v)
    a' = (a_inf - a)/tau_a
    b' = (b_inf - b)/tau_b
}

PROCEDURE rates(v (mV)) { :callable from hoc
    LOCAL alpha, beta
    TABLE a_inf, tau_a, b_inf, tau_b  : DEPEND vhalfn, vhalfl, tlmax, tnmax
    FROM -20 TO 130 WITH 750
    
    alpha = alpha_a(v)
    beta = beta_a(v)
    tau_a = 1/(alpha + beta)     
    a_inf = alpha*tau_a
    
    alpha = alpha_b(v)
    beta = beta_b(v)
    tau_b = 1/(alpha + beta)     
    b_inf = alpha*tau_b

}

