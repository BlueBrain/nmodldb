: Delayed rectifier K+ channel
: from Durstewitz & Gabriel (2006), Cerebral Cortex

NEURON {
	SUFFIX Kdr
	USEION k READ ki, ko WRITE ik
	RANGE gKdrbar, gk
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gKdrbar= 0.0338 (mho/cm2) <0,1e9>
}

ASSIGNED {
	v  (mV)
	ik (mA/cm2)
	gk (mho/cm2)
	ek (mV)
	ki (mM)
	ko (mM)
}

STATE {
	n 
}


INITIAL {
	n = alf(v)/(alf(v)+bet(v))
}

BREAKPOINT {
	SOLVE states METHOD derivimplicit
	gk = gKdrbar*n*n*n*n
	ek = 25*log(ko/ki)
	ik = gk*(v-ek)
}

DERIVATIVE states {
	n' = (1-n)*alf(v)-n*bet(v)
}

UNITSOFF

FUNCTION alf(v(mV)) (/ms) { 
	LOCAL va 
	va=v-13
	if (fabs(va)<1e-04) {
	  alf = -0.018*(-25 - va*0.5)
	}
	else {           
	  alf = -0.018*va/(-1+exp(-va/25))
	}
}

FUNCTION bet(v(mV)) (/ms) { 
	LOCAL vb 
	vb=v-23
	if (fabs(vb)<1e-04) {
	  bet =  0.0054*(12 - vb*0.5)
	} 
	else {
	  bet = 0.0054*vb/(-1+exp(vb/12))
	}
}	

UNITSON	
