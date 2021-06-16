TITLE Cardiac time independent inward rectifier IK1  current
:   from Courtemanche et al Am J Physiol 1998 275:H301


NEURON {
	SUFFIX IK1
	USEION k READ ek WRITE ik
	RANGE gK1, ik
	GLOBAL dummy : prevent vectorization for use with CVODE

}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
        (mM) = (milli/liter)
	 

}

PARAMETER {
	 gK1=0.1567e-3 (S/cm2) <0,1e9>
	
}


ASSIGNED {
	v (mV)
    	ik (mA/cm2)
	ek (mV)
	dummy
	      
}

LOCAL k
BREAKPOINT {
	
	ik = gK1/(1 + exp(0.07*(v + 80)))*(v - ek)
}





             
