TITLE slow potassium current

COMMENT Equations from 
   Golomb D, Amitai Y (1997) Propagating neuronal discharges in
   neocortical slices: computational and experimental study. J Neurophys
   78: 1199-1211.

>< Gating kinetics are at 36 degC. 
ENDCOMMENT

NEURON {
        SUFFIX Ks
        USEION k READ ek WRITE ik
        RANGE g, ik
	GLOBAL tau
}

UNITS {
	(S)  = (siemens)
        (mA) = (milliamp)
        (mV) = (millivolt)
}

PARAMETER {
        g 	(S/cm2)
	tau= 75 (ms)
}

ASSIGNED {
        v 	(mV)
	ek 	(mV)
        ik 	(mA/cm2)
	mtau 	(ms)
	minf
}

STATE { m }

BREAKPOINT { 
        SOLVE states METHOD cnexp
	ik= g* m* (v- ek) 
}

DERIVATIVE states {
	rates()
	m'= (minf- m)/ mtau
}

INITIAL {
	rates()
	m= minf 
}

PROCEDURE rates() { UNITSOFF
	minf= 1/ (1+ exp(-(v+ 39)/ 5))
	mtau= tau
} UNITSON

