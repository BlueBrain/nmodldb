TITLE  LG cell
:

NEURON {
	SUFFIX lg
	USEION K WRITE iK VALENCE 1
	USEION Na WRITE iNa VALENCE 1
	NONSPECIFIC_CURRENT iL
	GLOBAL gL,gbarNa ,gbarK,eL,eK, eNa
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
	gbarNa = 0.0035 (mho/cm2)
	gbarK = 0.004 (mho/cm2)
	gL = 0.00000729 (mho/cm2)
	eNa = 45 (mV)
	eK = -80 (mV)
	eL = -60 (mV)
}

STATE { Kon Koff Nahon Nahoff
}

ASSIGNED {
	iNa (mA/cm2)
	iK (mA/cm2)
	iL (mA/cm2)
}


BREAKPOINT {
	SOLVE ssss METHOD sparse
	iNa 	= 	gbarNa*minf(v)*minf(v)*minf(v)*Nahoff*(v-eNa)
	iK 	= 	gbarK*Koff*Koff*Koff*Koff*(v-eK)
	iL 	= 	gL*(v - eL)
}

KINETIC ssss {
	~ Nahon <-> Nahoff (hinf(v)/tauh(v),(1-hinf(v))/tauh(v))
	~ Kon <-> Koff (kinf(v)/tauk(v),(1-kinf(v))/tauk(v))
}


INITIAL {
		Nahon = hinf(v)
		Nahoff = 1 - Nahon
		Kon = kinf(v)
		Koff = 1 - Kon
}



FUNCTION minf(v){
	minf = 1/(1+exp(-0.08*(v+21)))
}

FUNCTION hinf(v){
	hinf = 1/(1+exp(0.13*(v+33)))
}

FUNCTION kinf(v){
	kinf = 1/(1+exp(-0.045*(v+33)))
}




FUNCTION tauh(v)(ms){
tauh =  (5/(1+exp(-0.12*(v+62))))
}




FUNCTION tauk(v(mV))(ms){
tauk =  (4 + 100/(1+exp(0.065*(v+5))))
}
