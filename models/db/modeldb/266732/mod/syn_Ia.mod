COMMENT
Tonic afferent inputs from muscle spindle
ENDCOMMENT
					       
NEURON {
	SUFFIX IaSyn
	RANGE gmax, e, i 
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
    gmax=0 (S/cm2)
	e=0	(mV)
}

ASSIGNED {
	v (mV)
	i (mA/cm2)
}

BREAKPOINT {
	i = gmax*(v - e)
}