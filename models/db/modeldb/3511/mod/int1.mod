TITLE  Int1 cell
:

NEURON {
	SUFFIX int1
	USEION K WRITE iK VALENCE 1
	USEION Na WRITE iNa VALENCE 1
	USEION h WRITE ih VALENCE 1
	NONSPECIFIC_CURRENT iL
	GLOBAL gL,gbarNa ,gbarK,gbarh,eL,eK, eNa,eh
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
	gbarNa = 0.0035 (mho/cm2)
	gbarK = 0.006 (mho/cm2)
	gL = 0.00000729 (mho/cm2)
	eh = -20 (mV)
	eNa = 45 (mV)
	eK = -80 (mV)
	eL = -30 (mV)
}

STATE { Kon Koff Nahon Nahoff hon hoff
}

ASSIGNED {
	ih (mA/cm2)
	iNa (mA/cm2)
	iK (mA/cm2)
	iL (mA/cm2)
}


BREAKPOINT {
	SOLVE ssss METHOD sparse
	ih = gbarh*hoff*(v-eh)
	iNa 	= 	gbarNa*minf(v)*minf(v)*minf(v)*Nahoff*(v-eNa)
	iK 	= 	gbarK*Koff*Koff*Koff*Koff*(v-eK)
	iL 	= 	gL*(v - eL)
}

KINETIC ssss {
	~ hon <-> hoff (hminf(v)/tauhm(v),(1-hminf(v))/tauhm(v))
	~ Nahon <-> Nahoff (hinf(v)/tauh(v),(1-hinf(v))/tauh(v))
	~ Kon <-> Koff (kinf(v)/tauk(v),(1-kinf(v))/tauk(v))
}


INITIAL {
                hon = hminf(v)
		hoff = 1 - hminf(v)
		Nahon = hinf(v)
		Nahoff = 1 - Nahon
		Kon = kinf(v)
		Koff = 1 - Kon
}

FUNCTION hminf(v){
	hminf = 1/(1+exp(2*(v+65)))
}


FUNCTION minf(v){
	minf = 1/(1+exp(-0.08*(v+26)))
}

FUNCTION hinf(v){
	hinf = 1/(1+exp(0.13*(v+38)))
}

FUNCTION kinf(v){
	kinf = 1/(1+exp(-0.045*(v+25)))
}




FUNCTION tauhm(v)(ms){
	tauhm = 200+2500/(1+exp(2*(v+65)))
}

FUNCTION tauh(v)(ms){
	tauh =  5/(1+exp(-0.12*(v+67)))
}




FUNCTION tauk(v(mV))(ms){
	tauk =  (4 + 150/(1+exp(0.065*(v+30))))
}
