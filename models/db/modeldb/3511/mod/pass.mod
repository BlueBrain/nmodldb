TITLE  Passive cell
:

NEURON {
	SUFFIX pass
	NONSPECIFIC_CURRENT iL
	RANGE gL,eL

}

UNITS {
	(mA) = (milliamp)
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(um) = (micron)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	v   (mV)
	dt  = 0.1 (ms)
	gL = 0.00000729 (mho/cm2)
	eL = -60 (mV)
}
STATE {
}
ASSIGNED {
	iL (mA/cm2)
}


BREAKPOINT {
	iL 	= 	gL*(v - eL)
}
