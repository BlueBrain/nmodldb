TITLE Transient Outward K+ current
COMMENT
	modified From DiFrancesco & Noble 1985 Phil Trans R Soc Lond 307:353-398 
    modified for Neuron by FE GANNIER
	francois.gannier@univ-tours.fr (University of TOURS)
ENDCOMMENT
INCLUDE "Unit.inc"
INCLUDE "Volume.inc"
NEURON {
	SUFFIX ito
	USEION k READ ek, ko, ki WRITE ik
	USEION ca READ cai
	RANGE g, ik, minf, mtau
}

PARAMETER {
	Kmca = 0.0005	(mM)
	Kmto = 10	(mM)
	g = 0.28	(uS/mM) 
:	kb = 4 (mM)
}

STATE { : r
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
	cai (mM)
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
:	ik = m * g/S * ((0.2 + ko) / (Km1 + ko)) * (cai / (Kmto + cai)) * ((v + 10) / (1 - exp(-0.2(/mV) * (v + 10)))) * (ki*exp(0.02(/mV) * v) - kb*exp(-0.02(/mV) * v))
:	correction
	ik = (1e-06)* m * g/S * (0.2 + (ko / (Kmto + ko))) * (cai / (Kmca + cai)) * (v + 10) / (1 - exp(-0.2(/mV) * (v + 10))) * (ki*exp(0.5 * v *F/RT) - ko*exp(-0.5 * v*F/RT))
}

DERIVATIVE states {
	rate(v)
	m' = (minf - m)/mtau
}

FUNCTION alp(v(mV)) (/ms) { 
	alp = (0.001)* 0.033(/s)*exp(-v / 17(mV))
}

FUNCTION bet(v(mV)) (/ms) { 
	bet = (0.001)* 33(/s) / (1 + exp(-(v + 10 (mV))/8(mV)))
}

: UNITSOFF
PROCEDURE rate(v(mV)) { LOCAL a, b
TABLE minf, mtau FROM -100 TO 100 WITH 200
	a = alp(v)  b = bet(v) 
	mtau = 1/(a + b)
	minf = a * mtau
}
: UNITSON 
