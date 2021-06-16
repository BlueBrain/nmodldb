TITLE pcell.mod   squid sodium, potassium, and leak channels
 
COMMENT
 Initialize this mechanism to steady-state voltage by calling
  rates_gsquid(v) from HOC, then setting m_gsquid=minf_gsquid, etc.
 Remember to set celsius=6.3 (or whatever) in your HOC file.
 See hh1.hoc for an example of a simulation using this model.
 SW Jaslove  6 March, 1992
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX pcell
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, gl, el
        GLOBAL minf, hinf, ninf, mexp, hexp, nexp
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 20 (degC)
        dt (ms)
        gnabar = .35 (mho/cm2)
        ena = 60 (mV)
        gkbar = .006 (mho/cm2)
        ek = -68 (mV)
        gl = .0005 (mho/cm2)
        el = -49 (mV)
}
 
STATE {
        m h n c
}
 
ASSIGNED {
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        minf hinf ninf mexp hexp nexp 
}
 
BREAKPOINT {
        SOLVE states
        ina = gnabar*m*m*m*m*h*(v - ena)
        ik = gkbar*n*n*(v - ek)      
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
 
PROCEDURE rates(v) {:Computes rate and o
         : ther constants at current v.
         : Call once from HOC to 
         : initialize inf at resting v.
     LOCAL  q10, tinc, alpha, beta, sum
     TABLE minf, mexp, hinf, hexp, ninf, nexp             DEPEND dt, celsius
FROM -100 TO 100 WITH 200
        q10 = 2.3^((celsius - 20)/10)
        tinc = -dt * q10
                :"m" sodium activation system
        alpha = .03 * vtrap(-(v+28),15)
        beta =  2.7 * exp(-(v+53)/18)
        sum = alpha + beta
        minf = alpha/sum
        mexp = 1 - exp(tinc*sum)
                :"h" sodium inactivation system
        alpha = .045 * exp(-(v+58)/18)
        beta = 0.72 / (exp(-(v+23)/14) + 1)
        sum = alpha + beta
        hinf = alpha/sum
        hexp = 1 - exp(tinc*sum)
                :"n" potassium activation system
        alpha = .024*vtrap(-(v-17),8) 
        beta = 0.2*exp(-(v+48)/35)
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

 