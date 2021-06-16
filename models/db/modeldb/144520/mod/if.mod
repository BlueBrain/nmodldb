TITLE Sodium/Calcium Cardiac hyperpolarizing IF current (HCN4)
COMMENT
	modified From DiFrancesco & Noble 1985 Phil Trans R Soc Lond 307:353-398 
    modified for Neuron by FE GANNIER
	francois.gannier@univ-tours.fr (University of TOURS)
ENDCOMMENT
INCLUDE "Unit.inc"
INCLUDE "Volume.inc"
NEURON {
	SUFFIX ifun
	USEION k READ ek, ko WRITE ik
	USEION na READ ena WRITE ina
	
	RANGE gK, gNa, ifun, ik, ina
	GLOBAL minf, mtau 
}

PARAMETER {
	gK = 3	(uS)
	gNa= 3	(uS)
	Kmf = 45	(mM)
}

STATE { : y
	m 
}

ASSIGNED {
	v (mV)
	celsius (degC) : 37
	ik (mA/cm2)
	ina (mA/cm2)
	ifun (mA/cm2)
	minf 
	mtau (ms)  
	ek (mV)
	ena (mV)
	ko (mM)
}

INITIAL {
	rate(v)
	m = minf
}

BREAKPOINT { LOCAL kc, ifk, ifna
	SOLVE states METHOD derivimplicit
: if = m * ifmax = m * (inafmax + ikfmax)
	kc = m * (ko/(ko+Kmf))
	ifna = (1e-06)*kc * (gNa/S*(v-ena))
	ifk  = (1e-06)*kc * (gK/S*(v-ek))
	ik = ifk
	ina = ifna
	ifun = ifk + ifna
}

DERIVATIVE states {
	rate(v)
	m' = (minf - m)/mtau
}

FUNCTION alp(v(mV)) (/ms) { 
: original
:	alp = (0.001)* 0.025(/s)*exp(-0.067(/mV)*(v + 52(mV)))
: correction
	alp = (0.001)* 0.05(/s)*exp(-0.067(/mV)*(v + 52(mV)-10))
}

FUNCTION bet(v(mV)) (/ms) { 
LOCAL Eo 
	Eo= v + 52 - 10
	if (fabs(Eo*1(/mV)) < 1e-5)
	{
		bet = (0.001)* 2.5 (/s)
	} else {
: original
:		bet = (0.001)* 0.5(/mV/s)*(v + 52(mV))/(1-exp(0.2(/mV)*(v + 52(mV))))
: correction Rq Eo = V + 42
		bet = (0.001)*1(/mV/s)*Eo/(1 - exp(-0.2(/mV)*Eo))
	}
}

: UNITSOFF
PROCEDURE rate(v(mV)) { LOCAL a,b :y
TABLE minf, mtau FROM -100 TO 100 WITH 200
	a = alp(v)  b = bet(v) 
	mtau = 1/(a + b)
	minf = a * mtau
}
: UNITSON 
