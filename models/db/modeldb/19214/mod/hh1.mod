TITLE hh1.mod   squid sodium, potassium, and leak channels

COMMENT
A replica of NEURON hh.mod membrane mechanism with 10 times 
smaller conductivities gnabar, gkbar and gl and slightly modified ek and el
for simulation of active dendritic membrane of the models described in: 
Korogod SM and Kulagina IB (1998) Biol Cybern 79:231-240
ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}

NEURON {
        SUFFIX hh1
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gna, gkbar, gk, gl, el
        GLOBAL minf, hinf, ninf, mexp, hexp, nexp
}

PARAMETER {
        v (mV)
        celsius = 6.3 (degC)
        dt (ms)
        gnabar = .012 (mho/cm2)
        ena = 50 (mV)
        gkbar = .0036 (mho/cm2)
        ek = -77.5 (mV)
        gl = .00003 (mho/cm2)
        el = -53.79 (mV)
}

STATE {
        m h n
}

ASSIGNED {
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        gna (mho/cm2)
        gk (mho/cm2)
        minf hinf ninf mexp hexp nexp
}

BREAKPOINT {
        SOLVE states
        gna = gnabar*m*m*m*h
        ina = gna*(v - ena)
        gk = gkbar*n*n*n*n
        ik = gk*(v - ek)
        il = gl*(v - el)
}

UNITSOFF

INITIAL {
   rates(v)
   m = minf
   h = hinf
   n = ninf
}

PROCEDURE states() {  :Computes state variables m, h, and n
        rates(v)      :             at the current v and dt.
        m = m + mexp*(minf-m)
        h = h + hexp*(hinf-h)
        n = n + nexp*(ninf-n)
}

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  q10, tinc, alpha, beta, sum
        TABLE minf, mexp, hinf, hexp, ninf, nexp DEPEND dt, celsius FROM -100 TO 100 WITH 200
        q10 = 3^((celsius - 6.3)/10)
        tinc = -dt * q10
                :"m" sodium activation system
        alpha = .1 * vtrap(-(v+40),10)
        beta =  4 * exp(-(v+65)/18)
        sum = alpha + beta
        minf = alpha/sum
        mexp = 1 - exp(tinc*sum)
                :"h" sodium inactivation system
        alpha = .07 * exp(-(v+65)/20)
        beta = 1 / (exp(-(v+35)/10) + 1)
        sum = alpha + beta
        hinf = alpha/sum
        hexp = 1 - exp(tinc*sum)
                :"n" potassium activation system
        alpha = .01*vtrap(-(v+55),10)
        beta = .125*exp(-(v+65)/80)
        sum = alpha + beta
        ninf = alpha/sum
        nexp = 1 - exp(tinc*sum)
}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}

UNITSON
