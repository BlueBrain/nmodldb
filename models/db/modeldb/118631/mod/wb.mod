TITLE Wang Buzsaki Model
:
:
:

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX wb
    USEION na READ ena WRITE ina
    USEION k READ ek WRITE ik
    RANGE gnabar, gkbar, vth
    RANGE m_inf, h_inf, n_inf
    RANGE tau_h, tau_n
    RANGE h_exp, n_exp
}


UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}

PARAMETER {
    gnabar  = .07  (mho/cm2)
    gkbar   = .09  (mho/cm2)

    ena = 55    (mV)
    ek  = -90   (mV)
    celsius = 36    (degC)
    dt              (ms)
    v               (mV)
    
    vth = 0
}

STATE {
    h n
}

ASSIGNED {
    ina (mA/cm2)
    ik  (mA/cm2)
    il  (mA/cm2)
    m_inf
    h_inf
    n_inf
    tau_h
    tau_n
    h_exp
    n_exp
    tadj
}


BREAKPOINT {
    SOLVE states
    ina = gnabar * m_inf*m_inf*m_inf*h * (v - ena)
    ik  = gkbar * n*n*n*n * (v - ek)
}


:DERIVATIVE states {   : exact Hodgkin-Huxley equations
:   evaluate_fct(v)
:   m' = (m_inf - m) / tau_m
:   h' = (h_inf - h) / tau_h
:   n' = (n_inf - n) / tau_n
:}

PROCEDURE states() {    : exact when v held constant
    evaluate_fct(v)
    h = h + h_exp * (h_inf - h)
    n = n + n_exp * (n_inf - n)
    VERBATIM
    return 0;
    ENDVERBATIM
}

UNITSOFF
INITIAL {
    h = 0
    n = 0
:
:  Q10 was assumed to be 3 for both currents
:
: original measurements at roomtemperature?
:    tadj = 3.0 ^ ((celsius-36)/ 10 )
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b
    tadj = 1.
    
    if (v == -35) { a = 1 
    }else{          a = -0.1 * (v+35.) / ( exp(-(v+35)/10.) - 1) }
    b = 4. * exp(-(v+60.)/18.)    
    m_inf = a / (a + b)

    a = 0.007 * exp(-0.05*(v+58.))
    b = 1. / ( 1 + exp(-(v+28)/10.) )
    tau_h = (1. / (a + b)) / tadj
    h_inf = a / (a + b)

    if (v == -34) { a = 0.1 
    }else{          a = -0.01 * (v+34.) / ( exp(-(v+34)/10.) - 1) }
    b = 0.125 * exp(-0.0125*(v+44))
    tau_n = (1. / (a + b)) / tadj
    n_inf = a / (a + b)

    h_exp = 1 - exp(-dt/tau_h)
    n_exp = 1 - exp(-dt/tau_n)
}

UNITSON
