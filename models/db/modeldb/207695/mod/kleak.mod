
COMMENT
K passive leak channel

ENDCOMMENT


NEURON {
    
	SUFFIX kleak
	
	USEION k READ ek WRITE ik
	RANGE g, ik, ek
}

PARAMETER {
	g = 0   	(S/cm2)
	

}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ik 	(mA/cm2)
    ek     (mV)
	v	(mV)
}
 

BREAKPOINT {

	ik = g * (v - ek)
} 




