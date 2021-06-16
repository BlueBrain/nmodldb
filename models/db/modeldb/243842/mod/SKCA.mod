
TITLE Calcium activated Potassium channel (SK)  

: Author: Chitaranjan Mahapatra (chitaranjan@iitb.ac.in)
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 

: For details refer: 
: Mahapatra C, Brain KL, Manchanda R, A biophysically constrained computational model of the action potential 
: of mouse urinary bladder smooth muscle. PLOS One (2018) 


:SK channel ( function of Calcium only)

NEURON {
	SUFFIX SKCA 
	USEION k READ ek WRITE ik
        USEION ca READ cai
        RANGE  gbar,gkahp,ik, ninf,taun,g
        GLOBAL Cq10,mt
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(pS) = (picosiemens)
	(um) = (micron)
}

PARAMETER {
	gbar = 0.01	(S/cm2)
       : n = 4
        cai = 50.e-6	(mM)
        a0 = 1.3e4	(1/ms-mM-mM-mM-mM)	
        b0 = 0.06	(1/ms)				
	    celsius = 37(degC)
	Cq10 = 2
	mt = 1
	
}

STATE {	n }

ASSIGNED {
	ik	(mA/cm2)
        g	(S/cm2)
        ninf
        taun	(ms)
	a	(1/ms)
        v	(mV)
        ek	(mV)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = gbar*n
	ik = g*(v-ek)
}

INITIAL {
	rate(cai)
	n=ninf
}

DERIVATIVE state {
	rate(cai)
	n' = (ninf - n)/taun
}


PROCEDURE rate(cai (mM)) {
	LOCAL q10
	q10 = Cq10^((celsius - 22 (degC))/10 (degC) )
	a = a0*cai^4
	taun = mt * q10* (0.2*(q10/(a + b0)))
	ninf = a/(a + b0)
}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	       trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}	