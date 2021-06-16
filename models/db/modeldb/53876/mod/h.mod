COMMENT
This file, h.mod, implements the hyperpol. activated (g_H) current from 
Quadroni and Knopfel 1994 table 1
ENDCOMMENT

NEURON {
	SUFFIX h
	NONSPECIFIC_CURRENT i
	RANGE i, Erev, gbar
	GLOBAL taun_fixed
}

UNITS {
	(S)	=	(siemens)
	(mV)	=	(millivolt)
	(mA)	=	(milliamp)
}

PARAMETER {
	gbar = 110e-6	(S/cm2) < 0, 1e9 >
	Erev = -46 (mV)
	taun_fixed = 50 (ms)
}

ASSIGNED {
	i (mA/cm2)
	v (mV)
	g (S/cm2)
	ninf
}

STATE {	n }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * n
	i = g * (v - Erev)
}

INITIAL {
	: assume that v has been constant for a long time
	n = alphan(v)/(alphan(v) + betan(v))
}

DERIVATIVE states {
	rates(v)
	n' = (ninf - n)/taun_fixed
}

FUNCTION alphan(Vm (mV)) (/ms) {
	UNITSOFF
	alphan = 0.02/(1 + exp((Vm + 90)/7.5))
	UNITSON
}

FUNCTION betan(Vm (mV)) (/ms) {
	UNITSOFF
	betan =  0.02/(1 + exp(-(Vm + 90)/7.5))
	UNITSON
}

PROCEDURE rates(Vm(mV)) {
	ninf = alphan(Vm)/(alphan(Vm) + betan(Vm))
}
