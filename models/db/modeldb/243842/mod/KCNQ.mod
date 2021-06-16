
TITLE Voltage dependent Potassium channel (KCNQ) 
 
: Author: Chitaranjan Mahapatra (chitaranjan@iitb.ac.in)
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 

: For details refer: 
: Mahapatra C, Brain KL, Manchanda R, A biophysically constrained computational model of the action potential 
: of mouse urinary bladder smooth muscle. PLOS One (2018) 

NEURON {
	SUFFIX KCNQ
	USEION k READ ek WRITE ik
        RANGE  gbar,ik,ap,vhalfl,k1,mt,ht,vhalf2,k2,ah
      GLOBAL hinf, htau, iinf,itau, q10
}
UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v 		(mV)
	ek
	gbar=0.009 	(mho/cm2)
    vhalfl=7.1   	(mV)
	kl=-15
    vhalf2=-54.5  	(mV)
	k2=9
    celsius (degC)
	q10 = 2.3
    ap = 2
	ah = 1
	mt = 1
	ht = 10
}

STATE {
         h i
}

ASSIGNED {
	ik (mA/cm2)
        iinf
	    itau
		hinf
		htau
       
}

INITIAL {
	rate(v)
	h= hinf
	i = iinf
	iinf = 0.5
}


BREAKPOINT {
	SOLVE state METHOD cnexp
	ik = gbar*h^ap*i^ah*(v-ek)
}

DERIVATIVE state {
        rate(v)
        i' = (iinf - i)/itau
		h' = (hinf-h)/htau
}

PROCEDURE rate(v (mV)) { :callable from hoc
        LOCAL a,qt, b
        qt=q10^((celsius-37)/10)
        hinf = (1/(1 + exp((v-vhalfl)/kl)))
		iinf = (0.55/(1 + exp((v-vhalf2)/k2)))+0.45
        htau = mt * (1/(1+(((v+15)/20)^2)))/qt
	    itau = ht * (200+ (10/(1+(((v+54.18)/120)^5))))/qt
}
FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	       trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}	