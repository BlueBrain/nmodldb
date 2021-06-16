COMMENT

Deterministic model of kinetics and voltage-dependence of Ih-currents
in layer 5 pyramidal neuron, see Kole et al., 2006. Implemented by
Stefan Hallermann.

Predominantly HCN1 / HCN2 

ENDCOMMENT

TITLE Ih-current

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
     (mM) = (milli/liter)

}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	dt 	   		(ms)
	v 	   		(mV)
        ehd=-47 		(mV) 				       
	ghdbar=0.00015 (mho/cm2)	
	gamma_ih	:not used
	seed		:not used
}


NEURON {
	SUFFIX ih
	NONSPECIFIC_CURRENT Iqq
	RANGE Iqq,ghdbar
}

STATE {
	qq
}

ASSIGNED {
	Iqq (mA/cm2)
}

INITIAL {
	qq=alpha(v)/(beta(v)+alpha(v))
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	Iqq = ghdbar*qq*(v-ehd)
}

FUNCTION alpha(v(mV)) {

	alpha = 0.001*6.43*(v+154.9)/(exp((v+154.9)/11.9)-1)
	: parameters are estimated by direct fitting of HH model to
        : activation time constants and voltage activation curve
        : recorded at 34C

}

FUNCTION beta(v(mV)) {
	beta = 0.001*193*exp(v/33.1)			
}

DERIVATIVE state {     : exact when v held constant; integrates over dt step
	qq' = (1-qq)*alpha(v) - qq*beta(v)
}
