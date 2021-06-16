TITLE Time-dependent (delayed) K+ current
COMMENT
	modified From DiFrancesco & Noble 1985 Phil Trans R Soc Lond 307:353-398 
    modified for Neuron by FE GANNIER
	francois.gannier@univ-tours.fr (University of TOURS)
ENDCOMMENT
INCLUDE "Unit.inc"
INCLUDE "Volume.inc"
NEURON {
	SUFFIX ikd
	USEION k READ ek, ki, ko WRITE ik
	RANGE ik, imax
	GLOBAL minf, mtau 
}

PARAMETER {
	imax = 180 (nA)
}

STATE { : x
	m 
}

ASSIGNED {
	v (mV)
	celsius (degC) : 37
	ik (mA/cm2)
	minf 
	mtau (ms)  
	ek (mV)
	ko (mM)
	ki (mM)
}

LOCAL RT
INITIAL {
	RT = (1000)*R*(273.15+celsius)
	rate(v)
	m = minf
}

BREAKPOINT { 
	SOLVE states METHOD derivimplicit
: original
:	ik = (1e-6) * m * imax/S * (ki - ko*exp(-v/25(mV)))/140(mM)
:	correction
	ik = (1e-06)* m * imax/S * (ki - ko*exp(-v*F/RT))/140(mM)
}

DERIVATIVE states {
	rate(v)
	m' = (minf - m)/mtau
}

FUNCTION alp(v(mV)) (/ms) { 
	alp = (0.001)* 0.5(/s)*exp(0.0826(/mV)*(v + 50(mV))) / (1 + exp(0.057(/mV)*(v + 50)))
}

FUNCTION bet(v(mV)) (/ms) { 
	bet = (0.001)* 1.3(/s)*exp(-0.06(/mV)*(v + 20(mV))) / (1 + exp(-0.04(/mV)*(v + 20 (mV))))
}

: UNITSOFF
PROCEDURE rate(v (mV)) { LOCAL a,b,c :
TABLE minf, mtau FROM -100 TO 100 WITH 200
	a = alp(v)  b = bet(v) 
	mtau = 1/(a + b)
	minf = a * mtau
}
: UNITSON 
