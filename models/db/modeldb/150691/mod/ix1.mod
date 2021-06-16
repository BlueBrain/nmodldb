TITLE Cardiac IKx1 current
: from BEELER & REUTER, J.Physiol, 1977

NEURON {
	SUFFIX IKx1
	USEION k READ ek WRITE ik
	RANGE gx1, ik, ix1bar, Tauact, minf, mtau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mM) = (milli/liter)
}

PARAMETER {
	gx1=0.0008 (S/cm2) <0,1e9>
	Tauact=1 (ms)
}

STATE { x1
	m 
}

ASSIGNED {
	v (mV)
	celsius (degC) : 37
	ik (mA/cm2)
	ix1bar (mA/cm2)
	minf 
	mtau (ms)  
	ek (mV)      
}

INITIAL {
	rate(v*1(/mV))
	m = minf
}

BREAKPOINT {
SOLVE states METHOD derivimplicit
	ix1bar = gx1*(exp(0.04*(v+ 77))-1)/exp(0.04*(v + 35))
	ik = ix1bar*m
}

DERIVATIVE states {	
	rate(v*1(/mV))
	m' = (minf - m)/mtau
}

UNITSOFF
FUNCTION alp(v(mV)) { 
	alp = 0.0005*exp(0.083*(v + 50))/(exp(0.057*(v + 50))+1)
}

FUNCTION bet(v(mV)) { 
	bet = 0.0013*exp(-0.06*(v + 20))/(exp(-0.04*(v + 20)) + 1)
}

PROCEDURE rate(v) 
{
LOCAL a,b
TABLE minf, mtau DEPEND celsius FROM -100 TO 100 WITH 200
	a = alp(v)  b = bet(v)
	mtau = 1/(a + b)
	minf = a/(a + b)
}
UNITSON