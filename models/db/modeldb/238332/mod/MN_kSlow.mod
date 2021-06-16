NEURON {
    SUFFIX MN_kSlow
    USEION k READ ek WRITE ik
    RANGE gmax, ik
    RANGE alpha_A, alpha_B, alpha_C, alpha_D, alpha_E
    RANGE beta_A, beta_B, beta_C, beta_D, beta_E
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
}

PARAMETER {
    gmax = 0  (S/cm2)

    alpha_A = 0.2
    alpha_B = 0.0
    alpha_C = 1.0
    alpha_D = -2.96
    alpha_E = -7.74
    
    beta_A = 0.05
    beta_B = 0.0
    beta_C = 1.0
    beta_D = -14.07
    beta_E = 6.1
}

ASSIGNED {
    v (mV)
    ek (mV)
    ik (mA/cm2)
    ks_ninf
    ks_ntau (ms)
}

STATE {
    n
}

INITIAL {
    ek = -80.0
    rates()
    n = ks_ninf
    
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik = gmax * n * (v - ek)
}

DERIVATIVE states {
    rates()
    n' =  ( ks_ninf  -  n ) /  ks_ntau
}

UNITSOFF
PROCEDURE rates() {LOCAL alpha_n, beta_n
    :alpha_n = 0.2/(1.0 + exp((-2.96+v)/-7.74))
    :beta_n = 0.05/(1.0 + exp((-14.07+v)/6.1))

    alpha_n = (alpha_A+alpha_B*v)/(alpha_C+exp((alpha_D+v)/alpha_E))
    beta_n = (beta_A+beta_B*v)/(beta_C+exp((beta_D+v)/beta_E))

    ks_ntau = 1/( alpha_n + beta_n)
    ks_ninf = alpha_n * ks_ntau
    
}
UNITSON


