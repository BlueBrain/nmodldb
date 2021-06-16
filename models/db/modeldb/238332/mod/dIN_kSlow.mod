NEURON {
    SUFFIX dIN_kSlow
    USEION k READ ek WRITE ik
    RANGE gmax, ik
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
}

PARAMETER {
    gmax = 0  (S/cm2)
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
    ek = -81.5
    rates()
    n = ks_ninf
    
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik = gmax * pow(n,2) * (v - ek)
}

DERIVATIVE states {
    rates()
    n' =  ( ks_ninf  -  n ) /  ks_ntau
}

UNITSOFF
PROCEDURE rates() {LOCAL alpha_n, beta_n
    alpha_n = alphabeta(4.61973318e-01,8.20458521e-03,4.59367292e+00,-4.20812882e+00,-1.19678988e+01)
    beta_n  = alphabeta(9.24268986e-02,-1.35395653e-03,1.61527685e+00,2.10656266e+05,3.32762273e+05)

    ks_ntau = 1/( alpha_n + beta_n)
    ks_ninf = alpha_n * ks_ntau
    
}

FUNCTION alphabeta(A,B,C,D,E){
    alphabeta = (A + B*v)/(C + exp((v+D)/E))
}
UNITSON


