TITLE  Int1s cell
:

NEURON {
	SUFFIX int1s
	USEION h WRITE ih VALENCE 1
	NONSPECIFIC_CURRENT iL
	GLOBAL gL,gbarh,eL,eh
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
	gbarh = 0.000 (mho/cm2)
	gL = 0.00000729 (mho/cm2)
	eh = -20 (mV)
	eL = -30 (mV)
}

STATE {hon hoff
}

ASSIGNED {
	ih (mA/cm2)
	iL (mA/cm2)
}


BREAKPOINT {
	SOLVE ssss METHOD sparse
	ih = gbarh*hoff*(v-eh)
	iL 	= 	gL*(v - eL)
}

KINETIC ssss {
	~ hon <-> hoff (hminf(v)/tauhm(v),(1-hminf(v))/tauhm(v))
}


INITIAL {
                hon = hminf(v)
}

FUNCTION hminf(v){
	hminf = 1/(1+exp(2*(v+55)))
}

FUNCTION tauhm(v)(ms){
	tauhm = 200+2500/(1+exp(2*(v+60)))
}
