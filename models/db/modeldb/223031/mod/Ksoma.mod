COMMENT

Potassium current for the soma
ENDCOMMENT
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX Ksoma
        USEION k READ ek WRITE ik
        RANGE gksoma, ik
        GLOBAL ninf, nexp, ntau
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius (degC)
        dt (ms)
        gksoma = .0319 (mho/cm2)
        ek (mV)
}
 
STATE {
        n 
}
 
ASSIGNED {
        ik (mA/cm2)
        ninf 
	nexp 
	ntau (ms)
}
 
INITIAL {
	rate(v)
	n = ninf
}

BREAKPOINT {
        SOLVE state METHOD cnexp
	ik = gksoma*n*n*n*n*(v - ek)    
}

DERIVATIVE state {	:exact when v held constant
	rate(v)
	n' = (ninf-n)/ntau
}
UNITSOFF
PROCEDURE rate(v(mV)) {  :Computes rate and other constants at 
		      :current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL q10, tinc, alpha, beta
        TABLE ninf, nexp, ntau DEPEND celsius FROM -200 TO 
100 WITH 300
		q10 = 3^((celsius - 24)/10)
		tinc = -dt*q10
		alpha = 0.018*vtrap(-(v-25),25)
		beta = 0.0036*vtrap(v-35,12)
		ntau = 1/(alpha + beta)
		ninf = alpha*ntau
		nexp = 1-exp(tinc/ntau)
}
FUNCTION vtrap(x,y) {	:Traps for 0 in denominator of rate eqns.
		if (fabs(x/y) < 1e-6) {
			vtrap = y*(1 - x/y/2)
		}else{
			vtrap = x/(exp(x/y) - 1)
		}
}
UNITSON
