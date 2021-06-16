NEURON {
    SUFFIX dIN_na
    USEION na READ ena WRITE ina
    RANGE gmax, ina
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
    ena (mV)
    ina (mA/cm2)
    na_minf
    na_hinf
    na_mtau (ms)
    na_htau (ms)
    rate_h
    rate_m
    
}

STATE {
    h 
    m
    
}

INITIAL {
    ena = 50.0
    rates()
    h = na_hinf
    m = na_minf
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ina = gmax * pow(m,3) * h * (v - ena)
}

DERIVATIVE states {
    rates()
    m' =  ( na_minf  -  m ) /  na_mtau
    h' =  ( na_hinf - h ) / na_htau
}

UNITSOFF
PROCEDURE rates() {LOCAL alpha_m, beta_m, alpha_h, beta_h
    alpha_m = alphabeta(8.67,0.0,1.0,-1.01,-12.56)
    beta_m =  alphabeta(3.82,0.0,1.0,9.01,9.69)
    alpha_h = alphabeta(0.08,0.0,0.0,38.88,26.0)
    beta_h =  alphabeta(4.08,0.0,1.0,-5.09,-10.21)

    na_mtau = 1/( alpha_m + beta_m)
    na_htau = 1/( alpha_h + beta_h)
    na_minf = alpha_m * na_mtau
    na_hinf = alpha_h * na_htau
}

FUNCTION alphabeta(A,B,C,D,E){
    alphabeta = (A + B*v)/(C + exp((v+D)/E))
}
UNITSON

