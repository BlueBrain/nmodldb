NEURON {
    SUFFIX dIN_kFast
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
    kf_ninf
    kf_ntau (ms)
     
}

STATE {
    n
}

INITIAL {
    ek = -81.5    
    rates()
    n = kf_ninf
    
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik = gmax * pow(n,4)  * (v - ek)  
}

DERIVATIVE states {
    rates()
    n' =  ( kf_ninf  -  n ) /  kf_ntau
}

UNITSOFF
PROCEDURE rates() {LOCAL alpha_n, beta_n
    alpha_n =alphabeta(5.05922619,0.0665406,5.12003207,-18.39568861,-25.42482239) 
    beta_n  =alphabeta(5.04899873e-01,0.0,0.0,2.86904073e+01,3.46245833e+01)

    kf_ntau = 1/( alpha_n + beta_n)
    kf_ninf = alpha_n * kf_ntau
}

FUNCTION alphabeta(A,B,C,D,E){
    alphabeta = (A + B*v)/(C + exp((v+D)/E))
}
UNITSON


