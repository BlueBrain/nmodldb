
COMMENT

gBoltzT.mod

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX Potassium
	USEION k READ ek WRITE ik
	RANGE n, gk, gbar, ninf, nexp
	GLOBAL v12, vSlope, tau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

PARAMETER {
	gbar = 150   	(pS/um2)	: 0.03 mho/cm2
	v 		(mV)
								
	v12 = -17.3	(mV)
	vSlope = 11.3	(mV)
	tau = 3
} 


ASSIGNED {
	ik 		(mA/cm2)
	gk		(pS/um2)
	ek		(mV)
	ninf
	nexp
}
 

STATE { n }

INITIAL { 
	states()
	
}

BREAKPOINT {
        SOLVE states
	gk = gbar*n
	ik = (1e-4) * gk * (v - ek)
} 

PROCEDURE states() {   : Computes state variable n 

	nexp = 1-exp(-dt/tau)
	ninf = 1/(1+exp(-(v-v12)/vSlope))
	n = n + nexp*(ninf-n)
}
