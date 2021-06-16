: Fast Na+ channel
: from Durstewitz & Gabriel (2006), Cerebral Cortex

NEURON {
	SUFFIX Naf
	USEION na WRITE ina
	RANGE gNafbar, gna, ena
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gNafbar= 0.086 (mho/cm2) <0,1e9>
        ena = 55 (mV)
}

ASSIGNED {
        v  (mV)
	ina (mA/cm2)
	gna (mho/cm2)
}

STATE {
	m 
        h 
}

INITIAL {
	m = malf(v)/(malf(v)+mbet(v))
	h = half(v)/(half(v)+hbet(v))
}

BREAKPOINT {
	SOLVE states METHOD derivimplicit
	gna = gNafbar*m*m*m*h
	ina = gna*(v-ena)
}

DERIVATIVE states {
	m' = (1-m)*malf(v)-m*mbet(v)
	h' = (1-h)*half(v)-h*hbet(v)
}

UNITSOFF

FUNCTION malf(v(mV)) (/ms) { 
	LOCAL va 
	va=v+28
	if (fabs(va)<1e-04){
	   malf = -0.2816*(-9.3-va*0.5)
	} 
	else {
	   malf = -0.2816*va/(-1+exp(-va/9.3))
	}
}

FUNCTION mbet(v(mV)) (/ms) { 
	LOCAL vb 
	vb=v+1
	if (fabs(vb)<1e-04) {
	   mbet = 0.2464*(6-vb*0.5)
	}
	else {
	   mbet = 0.2464*vb/(-1+exp(vb/6))
	}
}	

FUNCTION half(v(mV)) (/ms) { 
	half=0.098*exp(-(v+43.1)/20)
}

FUNCTION hbet(v(mV)) (/ms) { 
	hbet=1.4/(1+exp(-(v+13.1)/10)) 
}

UNITSON
