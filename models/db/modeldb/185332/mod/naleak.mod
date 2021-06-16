
COMMENT

Na passive leak channel

ENDCOMMENT


NEURON {
	SUFFIX naleak
	USEION na READ ena
	RANGE g, i, ena
    NONSPECIFIC_CURRENT i  
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
	i 	(mA/cm2)
    ena     (mV)
}
 

BREAKPOINT {

	i = g * (v - ena)
} 




