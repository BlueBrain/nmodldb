TITLE Persistent sodium current 

COMMENT Equations from 
   Golomb D, Amitai Y (1997) Propagating neuronal discharges in
   neocortical slices: computational and experimental study. J Neurophys
   78: 1199-1211.

>< Gating kinetics are at 36 degC. 
ENDCOMMENT

NEURON {
        SUFFIX NaP
        USEION na READ ena WRITE ina 
        RANGE g, ina
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
	ena	(mV)
        ina	(mA/cm2)
	minf
}

BREAKPOINT { 
	rates()
	ina= g* minf* (v- ena) 
}

INITIAL {
	rates()
}

PROCEDURE rates() { UNITSOFF
	minf= 1/ (1+ exp(-(v+ 40)/ 5))
} UNITSON

