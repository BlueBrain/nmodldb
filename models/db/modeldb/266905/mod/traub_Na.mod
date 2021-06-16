TITLE fast Na channel Traub

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
    (mS) = (millisiemens)
}

NEURON {
    THREADSAFE
    : note - The mechanism has the same name as the suffix
    : every variable accessible in NEURON will be having the suffix _gNa
    
    SUFFIX gNa
    USEION na READ ena WRITE ina
    RANGE gmax, g, i, tau_m, tau_h
    : GLOBAL m_inf, tau_m, h_inf, tau_h
}

PARAMETER {
    gmax=0.03 (mho/cm2)  : value at soma
}

STATE {
    m
    h
}

ASSIGNED {
    v 	(mV)
    ena (mV)
    ina		(mA/cm2)
    m_inf	(1)
    h_inf	(1)
    tau_m	(ms)
    tau_h	(ms)
    g		(S/cm2)
    i		(mA/cm2)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g = gmax*m^2*h
    i = g*(v-ena)
    ina = i
}

INITIAL {
    rates(v)
    m = m_inf
    h = h_inf
}


FUNCTION alpha_m(v(mV)) {
    alpha_m = 0.32*(13.1-v)/(exp( (13.1-v)/4 ) - 1)
}

FUNCTION beta_m(v(mV)) {
    beta_m = 0.28*(v - 40.1)/(exp( (v-40.1)/5 ) - 1)
}

FUNCTION alpha_h(v(mV)) {
    alpha_h = 0.128*exp( (17-v)/18 )
}

FUNCTION beta_h(v(mV)) {
    beta_h = 4/(1 + exp( (40 - v)/5 ) ) 
}

DERIVATIVE states {  
    rates(v)
    m' = (m_inf - m)/tau_m
    h' = (h_inf - h)/tau_h
}

PROCEDURE rates(v (mV)) { :callable from hoc
    LOCAL alpha, beta
    TABLE m_inf, tau_m, h_inf, tau_h  : DEPEND vhalfn, vhalfl, tlmax, tnmax
    FROM -20 TO 130 WITH 750
    
    alpha = alpha_m(v)
    beta = beta_m(v)
    tau_m = 1/(alpha + beta)     
    m_inf = alpha*tau_m
    
    alpha = alpha_h(v)
    beta = beta_h(v)
    tau_h = 1/(alpha + beta)     
    h_inf = alpha*tau_h

}

