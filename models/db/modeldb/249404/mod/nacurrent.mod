TITLE nacurrent.mod
 
COMMENT
The current implementation here is adapted from modeldb.yale.edu/136715 (the
Fleidervish et al 2010 author's hh_Cs_scaled.mod) which is a modified 
(wrt rate functions and temperature dependence) version of NEURON's hh.mod,
which is an implementation of the original Hodgkin-Huxley equations.

Note: This mechanism is temperature dependent.

Note: Unlike Hodgkin-Huxley (and hh.mod), this file only describes a
      sodium current.
ENDCOMMENT
 
UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
}
 
NEURON {
    SUFFIX nacurrent
    USEION na READ ena WRITE ina
    RANGE gnabar, gna
    GLOBAL minf, hinf, mtau, htau
    THREADSAFE
}
 
PARAMETER {
    gnabar = .0 (S/cm2)	<0, 1e9>
}
 
STATE {
    m (1)
    h (1)
}
 
ASSIGNED {
    v (mV)
    celsius (degC)
    ena (mV)
    gna (S/cm2)
    ina (mA/cm2)
    minf (1)
    hinf (1)
    mtau (ms)
    htau (ms)
}
 
BREAKPOINT {
    SOLVE states METHOD cnexp
    gna = gnabar * m * m * m * h
    ina = gna * (v - ena)
}
 
INITIAL {
    rates(v)
    m = minf
    h = hinf
}

DERIVATIVE states {  
    rates(v)
    m' = (minf - m) / mtau
    h' = (hinf - h) / htau
}
 
PROCEDURE rates(v(mV)) {
    : Computes rate and other constants at specified v.
    LOCAL  alpha, beta, sum, q10
    TABLE minf, mtau, hinf, htau DEPEND celsius FROM -100 TO 100 WITH 200

UNITSOFF
    q10 = 3 ^ ((celsius - 23) / 10)
    : "m" sodium activation 
    alpha = -.182 * vtrap(-(v + 40), 6)
    beta =  -.124 * vtrap((v + 40), 6)
    sum = alpha + beta
    mtau = 0.25 / (q10 * sum)
    minf = alpha / sum
    : "h" sodium inactivation
    alpha = -0.015 * vtrap((v + 66), 6)
    beta = -0.015 * vtrap(-(v + 66), 6)
    sum = alpha + beta
    htau = 1 / (q10 * sum)
    hinf = alpha / sum
}
 
FUNCTION vtrap(x, y) {
    : Avoids divide by zero errors in rate functions by replacing with limit
    if (fabs(x / y) < 1e-6) {
        vtrap = -y * (1 - x / y / 2)
    } else {
        vtrap = x / (1 - exp(x / y))
    }
}
 
UNITSON
