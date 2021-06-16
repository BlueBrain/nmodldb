TITLE Naf_No.mod   fast sodium channel
 
COMMENT
This is the original Hodgkin-Huxley treatment for the set of sodium channel found
in the squid giant axon membrane.
Some parameters have been changed to correspond to McIntyre and Grill (2002) "Extracellular
stimulation of central neurons"
Author: Balbi
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
		(S) = (siemens)
}
 
NEURON {
        SUFFIX Naf_No
        USEION na READ ena WRITE ina
        RANGE gnamax, gna
        RANGE minf, hinf, mtau, htau
}
 
PARAMETER { : value for Nodes
        gnamax = 0.2 (S/cm2)   <0,1e9>
}
 
STATE {
        m h
}
 
ASSIGNED {
        v (mV)
        ena (mV)

		gna (S/cm2)
        ina (mA/cm2)
        minf hinf
		mtau (ms) htau (ms)
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnamax*m*m*m*h
		ina = gna*(v - ena)
} 
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
}

DERIVATIVE states {
        rates(v)
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
}
 
PROCEDURE rates(v(mV)) {  
		:Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum

UNITSOFF
                :"m" sodium activation system
        alpha = 6.57 * vtrap(-10.4 - v, 10.3)
        beta =  .304 * vtrap(v + 15.7, 9.16)
        sum = alpha + beta
		mtau = 1/sum
        minf = alpha/sum
                :"h" sodium inactivation system
        alpha = 0.34 * vtrap(v+104, 11)
        beta =  12.6/(1+exp(-(v-21.8)/13.4))
        sum = alpha + beta
		htau = 1/sum
        hinf = alpha/sum
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON
