TITLE CA1NaP.mod   Persistent sodium current
 
COMMENT
Persistent Na for CA1 pyramid
Specifications in Lipowsky et al, J.Neurophys 76:2181-2191, 1996
Half-activation voltage can be varied (is -49mV in paper)
BPG 26-5-98
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX CA1NaP
        USEION na WRITE ina
        RANGE gnaP, ina, ena
        GLOBAL mPinf, mPexp
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        dt (ms)
        gnaP = 0.00017 (mho/cm2)
        ena = 65 (mV)
        vhalf = -49 (mV)
}
 
STATE {
        mP
}
 
ASSIGNED {
        ina (mA/cm2)
        mPinf mPexp
}
 
BREAKPOINT {
        SOLVE states
        ina = gnaP*mP*(v - ena)
}
 
UNITSOFF
 
INITIAL {
	rates(v)
	mP = mPinf
}

PROCEDURE states() {  :Computes state variable mP
        rates(v)      :             at the current v and dt.
        mP = mP + mPexp*(mPinf-mP)
}
 
PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize "inf" at resting v.
        LOCAL alpha, beta
        TABLE mPinf, mPexp DEPEND dt FROM -100 TO 50 WITH 150
                :"m" sodium activation system
        alpha = -1.74*(v-11)/(exp(-(v-11)/12.94)-1)
        beta = 0.06*(v-5.9)/(exp((v-5.9)/4.47)-1)
        mPinf = 1/(1+exp(-(v-vhalf)/5))
        mPexp = 1 - exp(-dt*(alpha+beta))
}
 
 
UNITSON

