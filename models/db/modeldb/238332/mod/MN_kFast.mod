NEURON {
    SUFFIX MN_kFast
    USEION k READ ek WRITE ik
    RANGE gmax, ik
    RANGE alpha_A, alpha_B, alpha_C, alpha_D, alpha_E
    RANGE beta_A, beta_B, beta_C, beta_D, beta_E
    RANGE kf_ninf, kf_ntau
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
}

PARAMETER {
    gmax = 0  (S/cm2)

    alpha_A = 3.1
    alpha_B = 0.0
    alpha_C = 1.0
    alpha_D = -27.5
    alpha_E = -9.3
    
    beta_A = 0.44
    beta_B = 0.0
    beta_C = 1.0
    beta_D = 8.98
    beta_E = 16.19
}

ASSIGNED {
    v (mV)
    ek (mV)
    ik (mA/cm2)
    kf_ninf
    kf_ntau (ms)
}

STATE {
    n
}

INITIAL {
    ek = -80.0
    rates()
    n = kf_ninf
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik = gmax * n  * (v - ek)  
}

DERIVATIVE states {
    rates()
    n' =  ( kf_ninf  -  n ) /  kf_ntau
}

UNITSOFF
PROCEDURE rates() {LOCAL alpha_n, beta_n
    :alpha_n = 3.1/(1.0 + exp((-27.5+v)/(-9.3)))
    :beta_n = 0.44/(1.0 + exp((8.98+v)/16.19))
    alpha_n = (alpha_A+alpha_B*v)/(alpha_C+exp((alpha_D+v)/alpha_E))
    beta_n = (beta_A+beta_B*v)/(beta_C+exp((beta_D+v)/beta_E))
    
    kf_ntau = 1/( alpha_n + beta_n)
    kf_ninf = alpha_n * kf_ntau
}
UNITSON


