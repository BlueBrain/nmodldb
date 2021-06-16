TITLE nadend.mod   squid sodium channels
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX nadend
        USEION na READ ena WRITE ina
        RANGE gnabar
        GLOBAL minf, hinf, mexp, hexp
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius (degC)
        dt (ms)
 	gnabar = 0.1250(mho/cm2) 
}
 
STATE {
        m h c
}
 
ASSIGNED {
	ena (mV) 
	ina (mA/cm2)
        minf hinf mexp hexp  
}
 
BREAKPOINT {
        SOLVE states
        ina = gnabar*m*m*m*m*h*(v - ena)
}
 
UNITSOFF
 
INITIAL {
     rates(v)
     m = minf
     h = hinf
}

PROCEDURE states() {  :Computes state variables m, h
        rates(v)      :             at the current v and dt.
        m = m + mexp*(minf-m)
        h = h + hexp*(hinf-h)
}
 
PROCEDURE rates(v) {:Computes rate and o
         : ther constants at current v.
         : Call once from HOC to 
         : initialize inf at resting v.
     LOCAL  q10, tinc, alpha, beta, sum
     TABLE minf, mexp, hinf, hexp   
	DEPEND dt, celsius

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
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON

 



