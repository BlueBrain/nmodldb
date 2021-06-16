TITLE A-type potassium current

COMMENT Equations from 
   Golomb D, Amitai Y (1997) Propagating neuronal discharges in
   neocortical slices: computational and experimental study. J Neurophys
   78: 1199-1211.

>< Gating kinetics are at 36 degC. 
ENDCOMMENT

NEURON {
        SUFFIX KA
        USEION k READ ek WRITE ik
        RANGE g, ik
}

UNITS {
	(S)  = (siemens)
        (mA) = (milliamp)
        (mV) = (millivolt)
}

PARAMETER {
        g	(S/cm2)
}

ASSIGNED {
        v	(mV)
	ek	(mV)
        ik	(mA/cm2)
	mtau	(ms)
	htau	(ms)
	minf
	hinf
}

STATE { h }

BREAKPOINT { 
        SOLVE states METHOD cnexp
	ik= g* minf^3* h* (v- ek) 
}

DERIVATIVE states {
	rates()
	h'= (hinf- h)/ htau
}

INITIAL {
	rates()
	h= hinf 
}

PROCEDURE rates() { UNITSOFF
	minf= 1/ (1+ exp(-(v+ 50)/ 20))
	hinf= 1/ (1+ exp((v+ 80)/ 6))
	htau= 15
} UNITSON

