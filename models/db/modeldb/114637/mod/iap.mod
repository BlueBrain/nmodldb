TITLE HH channel
: Mel-modified Hodgkin - Huxley conductances (after Ojvind et al.)
: Re-modified by Ojvind to exactly simulate the old version.  92-1-16
: (db) 22.12.97 modifications for CVode

NEURON {
	SUFFIX iapnew
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	GLOBAL inf
	RANGE gnabar, gkbar, ena, ek, ik
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
PARAMETER {
	v (mV)
	celsius = 37	(degC)
	gnabar=.20 (mho/cm2)
	gkbar=.12 (mho/cm2)
	ena = 40 (mV)
	ek = -85 (mV)
	naactvha = 40 (mV)
	naiactvha = 45 (mV)
	kactvha = 40 (mV)	
	naactslope = -3 (mV)
	nainactslope = 3 (mV)
	kactslope = -3 (mV)
	Nainactivationtau = 0.5
}

STATE {
	m h n
}

ASSIGNED {
	ina (mA/cm2)
	ik (mA/cm2)
	inf[3]
	tau[3]
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = gnabar*m*m*h*(v - ena)
	ik = gkbar*n*n*(v - ek)
}

DERIVATIVE states {	: exact when v held constant
	mhn(v*1(/mV))
	m' = (inf[0] - m)/tau[0]
	h' = (inf[1] - h)/tau[1]
	n' = (inf[2] - n)/tau[2]
}

FUNCTION varss(v, i) {
	if (i==0) {
		varss = 1 / (1 + exp((v + naactvha)/(naactslope))) :Na activation
	}
	else if (i==1) {
		varss = 1 / (1 + exp((v + naiactvha)/(nainactslope))) :Na inactivation
	}
	else {
		:varss = 0
		varss = 1 / (1 + exp((v + kactvha)/(kactslope))) :K activation
	}
}

FUNCTION vartau(i) {
	if (i==0) {
		vartau = 0.05  :Na activation tau
	}
	else if (i==1) {
		vartau = Nainactivationtau   :Na inactivation tau
	}
	else {
		vartau = 2     :K activation
	}
}

PROCEDURE mhn(v) {LOCAL a, b :rest = -70
	TABLE inf,tau 
	DEPEND celsius, dt, kactslope, kactvha, nainactslope,naiactvha,naactslope, naactvha, Nainactivationtau
	FROM -100 TO 100 WITH 2000 : .1 mV steps
	FROM i=0 TO 2 {
		tau[i] = vartau(i)
		inf[i] = varss(v,i)
	}
}
UNITSON
